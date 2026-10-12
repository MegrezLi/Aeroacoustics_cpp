#include "aeroacoustics.hpp"
#include "../src/acoustics/source_models.hpp"
#include "../src/turbine/analysis/arrival_history.hpp"
#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <new>
#include <numeric>

namespace {
bool counting = false;
std::size_t allocations = 0;
void check(bool condition, const char *message) {
    if (!condition) throw std::runtime_error(message);
}
void near(double a, double b, double tolerance = 1e-10) {
    check(a == b || (std::isfinite(a) && std::isfinite(b) && std::abs(a-b) <= tolerance),
          "Optimized value differs from reference");
}
}
void *operator new(std::size_t n) {
    auto p = std::malloc(n ? n : 1);
    if (!p) throw std::bad_alloc();
    if (counting) ++allocations;
    return p;
}
void *operator new[](std::size_t n) { return ::operator new(n); }
void operator delete(void *p) noexcept { std::free(p); }
void operator delete[](void *p) noexcept { std::free(p); }
void operator delete(void *p, std::size_t) noexcept { std::free(p); }
void operator delete[](void *p, std::size_t) noexcept { std::free(p); }

void storage() {
    using namespace aeroacoustics;
    counting = true;
    allocations = 0;
    detail::TnoWorkspace tno;
    counting = false;
    check(allocations == 0 && tno.factor.empty(), "Unused TNO allocated storage");
    tno.prepare_storage();
    counting = true;
    allocations = 0;
    tno.prepare_storage();
    counting = false;
    check(allocations == 0 && tno.factor.size() == 61*61, "TNO storage not reusable");

    Parameters p;
    std::vector<Node> nodes(6);
    for (std::size_t n=0; n<nodes.size(); ++n) nodes[n].aero_center = {0., double(n), 80.};
    std::vector<Vec3> observers(17);
    for (std::size_t o=0; o<observers.size(); ++o) observers[o] = {100., double(o), 2.};
    AcousticWorkspace full(p), blocked(p);
    const auto expected = full.evaluate(nodes, observers);
    ObserverBlockCallback consume = [&](std::size_t first, const Snapshot &block) {
        for (std::size_t i=0; i<block.size(); ++i)
            check(block[i] == expected[first+i], "Tail buffer changed output");
    };
    for (std::size_t block : {4, 1, 17, 20, 4}) {
        blocked.evaluate_blocks(nodes, observers, consume, block);
        counting = true;
        allocations = 0;
        for (int step=0; step<5; ++step) blocked.evaluate_blocks(nodes, observers, consume, block);
        counting = false;
        check(allocations == 0, "Repeated observer tail allocated memory");
        check(blocked.evaluate(nodes, observers) == expected, "Full snapshot after tail changed output");
    }
    auto copy = blocked;
    check(copy.evaluate(nodes, observers) == expected, "Workspace copy changed output");
}

void arrivals() {
    turbine::detail::ArrivalHistory packed;
    turbine::ArrivalSeries legacy;
    for (int i=0; i<5; ++i) {
        const std::vector<double> power{double(i*i), double(i+2)};
        packed.append(i*.7, i*.5, 8.+i, power);
        legacy.append({i*.7, i*.5, 8.+i, power});
    }
    auto independent = packed;
    packed.append(4., 3., 10., {10., 20.});
    std::size_t cursor = 0;
    std::size_t legacy_cursor = 0;
    turbine::ArrivalSample reused;
    for (double t : {0., .7, 1.01, 2.8, .2, 2.1}) {
        const auto expected = legacy.at(t);
        legacy.at_into(t, reused, legacy_cursor);
        near(reused.wind, 8.+t/.7);
        near(reused.emission, t*.5/.7);
        check(reused.power==expected.power, "Reusable arrival interpolation differs");
        std::vector<double> sum(2, 1.);
        independent.add_power(t, sum, cursor);
        for (std::size_t f=0; f<2; ++f) near(sum[f], 1.+expected.power[f], 0.);
    }
    std::vector<double> sum(2);
    cursor=0;
    independent.add_power(0., sum, cursor);
    allocations=0; counting=true;
    for (int i=0; i<1000; ++i) {
        std::fill(sum.begin(), sum.end(), 0.);
        independent.add_power(2.8*i/999, sum, cursor);
    }
    counting=false;
    check(allocations==0, "Packed interpolation allocated memory");
    allocations=0; counting=true;
    for(int i=0;i<1000;++i) legacy.at_into(2.8*i/999,reused,legacy_cursor);
    counting=false;
    check(allocations==0, "Reusable arrival interpolation allocated memory");
    bool rejected=false;
    try { independent.append(1., 4., 8., {0., 0.}); } catch (const std::invalid_argument &) { rejected=true; }
    check(rejected && independent.samples().size()==5, "Bad history append was accepted");
}

void turbulence() {
    using namespace aeroacoustics;
    double maximum=0.;
    // Include one-element windows, startup, repeated rebases and nearly constant flow.
    for (double dt : {5., .1, .00625}) for (int pattern=0; pattern<4; ++pattern) {
        TurbulenceState state({5.}, 1, dt, 100., 2);
        const auto samples=std::size_t(std::floor(5/dt+.5));
        Spectrum ring(samples), speed{50.};
        std::vector<Vec3> wind(1), leading{{0., 0., 100.}};
        for (std::size_t i=1; i<=samples*5+100; ++i) {
            const double v=pattern==0 ? 8. : pattern==1 ? 8.+1e-8*std::sin(i*.7) :
                           pattern==2 ? 10.+3.*std::sin(i*.013) :
                           (i<samples*2 ? 2.+double(i%7) : 5.);
            wind[0]={v,0.,0.};
            const double value=std::sqrt(v*v);
            ring[i<=samples ? i-1 : i%samples]=value;
            double expected=0.;
            if (i>samples) {
                const double mean=std::accumulate(ring.begin(),ring.end(),0.)/samples;
                double variance=0.;
                for(double x:ring) variance+=(x-mean)*(x-mean);
                expected=std::sqrt(variance/samples)/mean;
            }
            state.update(speed,wind,leading);
            maximum=std::max(maximum,std::abs(state.values[0]-expected));
            near(state.values[0],expected,1e-10);
            if(i==samples+2) {
                auto clone=state;
                clone.update(speed,wind,leading);
                auto copied_ring=ring;
                copied_ring[(i+1)%samples]=value;
                const double mean=std::accumulate(copied_ring.begin(),copied_ring.end(),0.)/samples;
                double variance=0.;
                for(double x:copied_ring) variance+=(x-mean)*(x-mean);
                near(clone.values[0],std::sqrt(variance/samples)/mean,1e-10);
                near(state.values[0],expected,1e-10);
            }
        }
    }
    std::cout << "TI maximum absolute difference: " << maximum << '\n';
}
int main() {
    try {
        storage(); arrivals(); turbulence();
        std::cout << "Optimization checks passed\n";
    } catch(const std::exception &e) {
        counting=false;
        std::cerr << e.what() << '\n';
        return 1;
    }
}
