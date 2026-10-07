#include "aeroacoustics_c.h"
#include "acoustic_levels.hpp"
#include "kernels.hpp"
#include <algorithm>
#include <atomic>
#include <exception>
#include <stdexcept>
namespace {
thread_local std::string error;
std::atomic<unsigned long long> calls{0};

void require(bool condition, const std::string &message) {
    if (!condition)
        throw std::invalid_argument(message);
}
void validate_inputs(int op, int n, const double *f, const double *x, const int *flags) {
    require(op >= 1 && op <= 7, "Unknown acoustic kernel operation");
    for (int i = 0; i < n; ++i)
        require(std::isfinite(f[i]) && f[i] > 0,
                "Frequency must be finite and positive at index " + std::to_string(i));
    auto finite = [&](int i) {
        require(std::isfinite(x[i]), "Non-finite C ABI x[" + std::to_string(i) + "]");
    };
    auto positive = [&](int i) {
        finite(i);
        require(x[i] > 0, "C ABI x[" + std::to_string(i) + "] must be positive");
    };
    auto nonnegative = [&](int i) {
        finite(i);
        require(x[i] >= 0, "C ABI x[" + std::to_string(i) + "] must be nonnegative");
    };
    positive(2); // Raw kernels divide by speed; unlike the driver they do not clamp it.
    if (op != 7)
        positive(1);
    if (op != 6) {
        positive(16);
        require(x[2] < x[16], "C ABI speed must be subsonic");
        finite(3);
        finite(4);
        positive(6);
        if (op != 3)
            positive(5);
        if (op != 7)
            finite(0); // Radians for op 4, degrees for the BPM kernels.
    }
    const bool bpm = op == 1 || op == 2 || op == 5;
    if (bpm) {
        require(flags[0] == 1 || flags[0] == 2, "Invalid BL method");
        if (flags[0] == 1)
            require(flags[1] >= 0 && flags[1] <= 2, "Invalid trip flag");
        if (flags[0] == 1 || op == 2)
            finite(7);
        if (flags[0] == 2) {
            if (op == 1)
                positive(8);
            else {
                positive(9);
                positive(10);
            }
        }
        if (flags[0] == 1 || op == 1 || op == 2)
            positive(17);
    }
    if (op == 3) {
        positive(20);
        require(flags[2] == 0 || flags[2] == 1, "Rounded-tip flag must be 0 or 1");
    }
    if (op == 4) {
        nonnegative(11);
        positive(18);
        positive(19);
    }
    if (op == 5) {
        positive(12);
        nonnegative(13);
        require(x[13] <= 14, "Bluntness angle must be within [0,14] degrees");
    }
    if (op == 6) {
        nonnegative(14);
        nonnegative(15);
    }
    if (op == 7) {
        positive(17);
        positive(18);
        // Nonpositive Cf disables that side in the legacy TNO kernel. Signed
        // edge ratios retain their magnitude convention; inactive sides may be zero.
        for (int side = 0; side < 2; ++side) {
            const int thickness = side == 0 ? 23 : 8;
            finite(thickness);
            finite(21 + side);
            finite(24 + side);
            if (x[21 + side] > 0) {
                positive(thickness);
                require(x[24 + side] != 0, "Active TNO side requires nonzero edge velocity ratio");
            }
        }
    }
}
} // namespace
extern "C" {
const char *aeroacoustics_last_error() { return error.c_str(); }
const char *aeroacoustics_backend() {
    static std::string name = aeroacoustics::backend();
    return name.c_str();
}
unsigned long long aeroacoustics_call_count() { return calls.load(); }
int aeroacoustics_kernel(int op, int n, const double *f, const double *x, const int *flags, double *out) {
    try {
        using namespace aeroacoustics;
        require(n > 0 && f && x && flags && out, "Invalid C ABI buffer");
        validate_inputs(op, n, f, x, flags);
        Parameters p;
        p.freqlist.assign(f, f + n);
        // Unused caller slots may contain sentinels. In particular, TNO's
        // Parameters validation must not inspect fields unused by operation 7.
        if (op != 6)
            p.spdsound = x[16];
        if (op == 1 || op == 2 || (op == 5 && flags[0] == 1) || op == 7)
            p.kinvisc = x[17];
        if (op == 4 || op == 7)
            p.airdens = x[18];
        if (op == 4)
            p.lturb = x[19];
        if (op == 1 || op == 2 || op == 5) {
            p.x_blmethod = flags[0];
            p.itrip = flags[1];
        }
        if (op == 3)
            p.round = flags[2] != 0;
        Spectrum a, b, c;
        switch (op) {
        case 1:
            a = lblvs(x[0], x[1], x[2], x[3], x[4], x[5], x[6], p, x[8], x[9], x[10], x[7]);
            break;
        case 2:
            std::tie(a, b, c) = tblte(x[0], x[1], x[2], x[3], x[4], x[5], x[6], p, x[8], x[9], x[10], x[7]);
            break;
        case 3:
            a = tipnois(x[0], x[20], x[1], x[2], x[3], x[4], x[6], p);
            break;
        case 4:
            a = inflownoise(x[0], x[1], x[2], x[3], x[4], x[5], x[6], x[11], p);
            break;
        case 5:
            a = blunt(x[0], x[1], x[2], x[3], x[4], x[5], x[6], x[12], x[13], p, x[8], x[9], x[10], x[7]);
            break;
        case 6:
            a = simple_guidati(x[2], x[1], x[15], x[14], p);
            break;
        case 7: {
            BoundaryLayer bl;
            bl.cf = {{x[21], x[22]}};
            bl.d99 = {{x[23], x[8]}};
            bl.edge_velocity_ratio = {{x[24], x[25]}};
            std::tie(a, b) = tblte_tno(x[2], x[3], x[4], x[5], x[6], bl, p);
            break;
        }
        }
        // Validate every channel before touching caller output, including op 2/7.
        for (const auto *v : {&a, &b, &c})
            for (double value : *v)
                validate_level(value);
        std::copy(a.begin(), a.end(), out);
        if (!b.empty())
            std::copy(b.begin(), b.end(), out + std::size_t(n));
        if (!c.empty())
            std::copy(c.begin(), c.end(), out + 2 * std::size_t(n));
        ++calls;
        error.clear();
        return 0;
    } catch (const std::exception &e) {
        error = e.what();
        return 1;
    } catch (...) {
        error = "Unknown C++ failure";
        return 2;
    }
}
}
