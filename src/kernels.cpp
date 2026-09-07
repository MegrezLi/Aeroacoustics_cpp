// Numerical expressions ported from OpenFAST (Apache-2.0), see NOTICE.
// Direct C++ formulas; validated against the original Fortran numerical body.
#include "kernels.hpp"
#include <algorithm>
#include <cmath>
namespace aeroacoustics {
constexpr double aa_epsilon=1e-16, pi=3.14159265358979323846, twopi=2*pi;
constexpr int x_blmethod_tables=2, itrip_none=0, itrip_heavy=1, itrip_light=2;
double log10aa(double x) {
    double f = 0;
    f = std::log10(std::max<double>(aa_epsilon, x));
    return f;
}

Spectrum lblvs(double alpstar, double c, double u, double theta, double phi, double l, double r, const Parameters& p, double d99var2, double dstarvar1, double dstarvar2, double stallval) {
    double d = 0;
    double dbarh = 0;
    double deltap = 0;
    double dstrp = 0;
    double dstrs = 0;
    double e = 0;
    double g1 = 0;
    double g2 = 0;
    double g3 = 0;
    int i = 0;
    double m = 0;
    double rc = 0;
    double rc0 = 0;
    double scale = 0;
    double st1prim = 0;
    double stpkprm = 0;
    double stprim = 0;
    Spectrum spllam(p.freqlist.size(), 0.);
    m = (u / static_cast<double>(p.spdsound));
    rc = ((u * c) / static_cast<double>(p.kinvisc));
    if ((p.x_blmethod == x_blmethod_tables)) {
        deltap = d99var2;
        dstrs = dstarvar1;
        dstrp = dstarvar2;
    } else {
        std::tie(deltap, dstrs, dstrp) = thick(c, rc, alpstar, p, stallval);
    }
    dbarh = directh_te(m, theta, phi);
    if ((dbarh <= 0)) {
        std::fill(spllam.begin(), spllam.end(), 0.0);
        return spllam;
    }
    if ((rc <= 130000.0)) {
        st1prim = 0.18;
    } else {
        if ((rc <= 400000.0)) {
            st1prim = (0.001756 * std::pow(rc, 0.3931));
        } else {
            st1prim = 0.28;
        }
    }
    stpkprm = (std::pow(10.0, (-(0.04) * alpstar)) * st1prim);
    if ((alpstar <= 3.0)) {
        rc0 = std::pow(10.0, ((0.215 * alpstar) + 4.978));
    } else {
        rc0 = std::pow(10.0, ((0.12 * alpstar) + 5.263));
    }
    d = (rc / static_cast<double>(rc0));
    if ((d <= 0.3237)) {
        g2 = ((77.852 * log10aa(d)) + 15.328);
    } else {
        if ((d <= 0.5689)) {
            g2 = ((65.188 * std::log10(d)) + 9.125);
        } else {
            if ((d <= 1.7579)) {
                g2 = (-(114.052) * std::pow(std::log10(d), 2));
            } else {
                if ((d <= 3.0889)) {
                    g2 = ((-(65.188) * std::log10(d)) + 9.125);
                } else {
                    g2 = ((-(77.852) * std::log10(d)) + 15.328);
                }
            }
        }
    }
    g3 = (171.04 - (3.03 * alpstar));
    scale = (10.0 * log10aa(((((deltap * std::pow(m, 5)) * dbarh) * l) / static_cast<double>(std::pow(r, 2)))));
    for (i = 1; i < (p.freqlist.size() + 1); ++i) {
        stprim = ((p.freqlist[(i - 1)] * deltap) / static_cast<double>(u));
        e = (stprim / static_cast<double>(stpkprm));
        if ((e <= 0.5974)) {
            g1 = ((39.8 * log10aa(e)) - 11.12);
        } else {
            if ((e <= 0.8545)) {
                g1 = ((98.409 * std::log10(e)) + 2.0);
            } else {
                if ((e <= 1.17)) {
                    g1 = (-(5.076) + std::sqrt((2.484 - (506.25 * std::pow(std::log10(e), 2)))));
                } else {
                    if ((e <= 1.674)) {
                        g1 = ((-(98.409) * std::log10(e)) + 2.0);
                    } else {
                        g1 = ((-(39.8) * std::log10(e)) - 11.12);
                    }
                }
            }
        }
        spllam[(i - 1)] = (((g1 + g2) + g3) + scale);
    }
    return spllam;
}

std::tuple<Spectrum, Spectrum, Spectrum> tblte(double alpstar, double c, double u, double theta, double phi, double l, double r, const Parameters& p, double d99var2, double dstarvar1, double dstarvar2, double stallval) {
    double a = 0;
    double a0 = 0;
    double a02 = 0;
    double aa = 0;
    double amaxa = 0;
    double amaxa0 = 0;
    double amaxa02 = 0;
    double amaxb = 0;
    double amina = 0;
    double amina0 = 0;
    double amina02 = 0;
    double aminb = 0;
    double ara0 = 0;
    double ara02 = 0;
    double b = 0;
    double b0 = 0;
    double bb = 0;
    double beta = 0;
    double beta0 = 0;
    double bmaxb = 0;
    double bmaxb0 = 0;
    double bminb = 0;
    double bminb0 = 0;
    double brb0 = 0;
    double dbarh = 0;
    double dbarl = 0;
    double delk1 = 0;
    double deltap = 0;
    double dstrp = 0;
    double dstrs = 0;
    double gamma = 0;
    double gamma0 = 0;
    int i = 0;
    double k1 = 0;
    double k2 = 0;
    double m = 0;
    double rc = 0;
    double rdstrp = 0;
    double rdstrs = 0;
    double separated = 0;
    double st1 = 0;
    double st1prim = 0;
    double st2 = 0;
    double stp = 0;
    double stpeak = 0;
    double sts = 0;
    double xcheck = 0;
    Spectrum splp(p.freqlist.size(), 0.);
    Spectrum spls(p.freqlist.size(), 0.);
    Spectrum splalph(p.freqlist.size(), 0.);
    m = (u / static_cast<double>(p.spdsound));
    rc = ((u * c) / static_cast<double>(p.kinvisc));
    if ((p.x_blmethod == x_blmethod_tables)) {
        deltap = d99var2;
        dstrs = dstarvar1;
        dstrp = dstarvar2;
    } else {
        std::tie(deltap, dstrs, dstrp) = thick(c, rc, alpstar, p, stallval);
    }
    dbarl = directl(m, theta, phi);
    dbarh = directh_te(m, theta, phi);
    rdstrs = ((dstrs * u) / static_cast<double>(p.kinvisc));
    rdstrp = ((dstrp * u) / static_cast<double>(p.kinvisc));
    st1 = (0.02 * std::pow(m, -(0.6)));
    if ((alpstar <= 1.333)) {
        st2 = st1;
    } else {
        if ((alpstar <= stallval)) {
            st2 = (st1 * std::pow(10.0, (0.0054 * std::pow((alpstar - 1.333), 2))));
        } else {
            st2 = (4.72 * st1);
        }
    }
    st1prim = ((st1 + st2) / static_cast<double>(2.0));
    a0 = a0comp(rc);
    a02 = a0comp((3.0 * rc));
    amina0 = amin(a0);
    amaxa0 = amax(a0);
    amina02 = amin(a02);
    amaxa02 = amax(a02);
    ara0 = ((20.0 + amina0) / static_cast<double>((amina0 - amaxa0)));
    ara02 = ((20.0 + amina02) / static_cast<double>((amina02 - amaxa02)));
    if ((rc < 95200.0)) {
        b0 = 0.3;
    } else {
        if ((rc < 857000.0)) {
            b0 = ((-(4.48e-13) * std::pow((rc - 857000.0), 2)) + 0.56);
        } else {
            b0 = 0.56;
        }
    }
    bminb0 = bmin(b0);
    bmaxb0 = bmax(b0);
    brb0 = ((20.0 + bminb0) / static_cast<double>((bminb0 - bmaxb0)));
    stpeak = st1;
    if ((rc < 247000.0)) {
        k1 = ((-(4.31) * log10aa(rc)) + 156.3);
    } else {
        if ((rc <= 800000.0)) {
            k1 = ((-(9.0) * std::log10(rc)) + 181.6);
        } else {
            k1 = 128.5;
        }
    }
    if ((rdstrp <= 5000.0)) {
        delk1 = (-(alpstar) * (5.29 - (1.43 * log10aa(rdstrp))));
    } else {
        delk1 = 0.0;
    }
    gamma = ((27.094 * m) + 3.31);
    beta = ((72.65 * m) + 10.74);
    gamma0 = ((23.43 * m) + 4.651);
    beta0 = ((-(34.19) * m) - 13.82);
    if ((alpstar <= (gamma0 - gamma))) {
        k2 = -(1000.0);
    } else {
        if ((alpstar <= (gamma0 + gamma))) {
            k2 = (std::sqrt((std::pow(beta, 2) - (std::pow((beta / static_cast<double>(gamma)), 2) * std::pow((alpstar - gamma0), 2)))) + beta0);
        } else {
            k2 = -(12.0);
        }
    }
    k2 = (k2 + k1);
    xcheck = gamma0;
    separated = false;
    if (((alpstar >= xcheck) || (alpstar > stallval))) {
        separated = true;
    }
    for (i = 1; i < (p.freqlist.size() + 1); ++i) {
        stp = ((p.freqlist[(i - 1)] * dstrp) / static_cast<double>(u));
        a = log10aa((stp / static_cast<double>(stpeak)));
        amina = amin(a);
        amaxa = amax(a);
        aa = (amina + (ara0 * (amaxa - amina)));
        splp[(i - 1)] = ((((aa + k1) - 3.0) + (10.0 * log10aa(((((dstrp * std::pow(m, 5)) * dbarh) * l) / static_cast<double>(std::pow(r, 2)))))) + delk1);
        sts = ((p.freqlist[(i - 1)] * dstrs) / static_cast<double>(u));
        if (!(separated)) {
            a = log10aa((sts / static_cast<double>(st1prim)));
            amina = amin(a);
            amaxa = amax(a);
            aa = (amina + (ara0 * (amaxa - amina)));
            spls[(i - 1)] = (((aa + k1) - 3.0) + (10.0 * log10aa(((((dstrs * std::pow(m, 5)) * dbarh) * l) / static_cast<double>(std::pow(r, 2))))));
            b = log10aa((sts / static_cast<double>(st2)));
            bminb = bmin(b);
            bmaxb = bmax(b);
            bb = (bminb + (brb0 * (bmaxb - bminb)));
            splalph[(i - 1)] = ((bb + k2) + (10.0 * log10aa(((((dstrs * std::pow(m, 5)) * dbarh) * l) / static_cast<double>(std::pow(r, 2))))));
        } else {
            spls[(i - 1)] = (10.0 * log10aa(((((dstrs * std::pow(m, 5)) * dbarl) * l) / static_cast<double>(std::pow(r, 2)))));
            splp[(i - 1)] = (10.0 * log10aa(((((dstrp * std::pow(m, 5)) * dbarl) * l) / static_cast<double>(std::pow(r, 2)))));
            b = log10aa((sts / static_cast<double>(st2)));
            aminb = amin(b);
            amaxb = amax(b);
            bb = (aminb + (ara02 * (amaxb - aminb)));
            splalph[(i - 1)] = ((bb + k2) + (10.0 * log10aa(((((dstrs * std::pow(m, 5)) * dbarl) * l) / static_cast<double>(std::pow(r, 2))))));
        }
        if ((splp[(i - 1)] < -(100.0))) {
            splp[(i - 1)] = -(100.0);
        }
        if ((spls[(i - 1)] < -(100.0))) {
            spls[(i - 1)] = -(100.0);
        }
        if ((splalph[(i - 1)] < -(100.0))) {
            splalph[(i - 1)] = -(100.0);
        }
    }
    return {splp, spls, splalph};
}

Spectrum tipnois(double alphtip, double alprat2, double c, double u, double theta, double phi, double r, const Parameters& p) {
    double alptipp = 0;
    double dbarh = 0;
    int i = 0;
    double l = 0;
    double m = 0;
    double mm = 0;
    double scale = 0;
    double stpp = 0;
    double term = 0;
    double um = 0;
    Spectrum spltip(p.freqlist.size(), 0.);
    if ((alphtip == 0.0)) {
        std::fill(spltip.begin(), spltip.end(), 0);
        return spltip;
    } else {
        if ((alphtip < 0.0)) {
            ;
        }
    }
    alptipp = (std::abs(alphtip) * alprat2);
    m = (u / static_cast<double>(p.spdsound));
    dbarh = directh_te(m, theta, phi);
    if (p.round) {
        l = ((0.008 * alptipp) * c);
    } else {
        if ((std::abs(alptipp) <= 2.0)) {
            l = ((0.023 + (0.0169 * alptipp)) * c);
        } else {
            l = ((0.0378 + (0.0095 * alptipp)) * c);
        }
    }
    mm = ((1.0 + (0.036 * alptipp)) * m);
    um = (mm * p.spdsound);
    term = (((((m * m) * std::pow(mm, 3)) * std::pow(l, 2)) * dbarh) / static_cast<double>(std::pow(r, 2)));
    if ((term != 0.0)) {
        scale = (10.0 * std::log10(term));
    } else {
        scale = 0.0;
    }
    for (i = 1; i < (p.freqlist.size() + 1); ++i) {
        stpp = ((p.freqlist[(i - 1)] * l) / static_cast<double>(um));
        spltip[(i - 1)] = ((126.0 - (30.5 * std::pow((log10aa(stpp) + 0.3), 2))) + scale);
    }
    return spltip;
}

Spectrum inflownoise(double alphanoise, double chord, double u, double theta, double phi, double d, double robs, double tinoise, const Parameters& p) {
    double beta2 = 0;
    double dbarh = 0;
    double dbarl = 0;
    double directivity = 0;
    double frequency_cutoff = 0;
    int i = 0;
    double kbar = 0;
    double ke = 0;
    double khat = 0;
    double lfc = 0;
    double mach = 0;
    double sears = 0;
    double splhigh = 0;
    double tinooisess = 0;
    double wavenumber = 0;
    Spectrum splti(p.freqlist.size(), 0.);
    mach = (u / static_cast<double>(p.spdsound));
    tinooisess = tinoise;
    dbarl = directl(mach, theta, phi);
    dbarh = directh_le(mach, theta, phi);
    if ((dbarh <= 0)) {
        std::fill(splti.begin(), splti.end(), 0.0);
        return splti;
    }
    frequency_cutoff = (((10 * u) / static_cast<double>(pi)) / static_cast<double>(chord));
    ke = (3.0 / static_cast<double>((4.0 * p.lturb)));
    beta2 = (1 - (mach * mach));
    for (i = 1; i < (p.freqlist.size() + 1); ++i) {
        if ((p.freqlist[(i - 1)] <= frequency_cutoff)) {
            directivity = dbarl;
        } else {
            directivity = dbarh;
        }
        wavenumber = ((twopi * p.freqlist[(i - 1)]) / static_cast<double>(u));
        kbar = ((wavenumber * chord) / static_cast<double>(2.0));
        khat = (wavenumber / static_cast<double>(ke));
        splhigh = ((10.0 * log10aa((((((((((std::pow(p.airdens, 2) * std::pow(p.spdsound, 4)) * p.lturb) * (d / static_cast<double>(2.0))) / static_cast<double>(std::pow(robs, 2))) * std::pow(mach, 5)) * std::pow(tinooisess, 2)) * std::pow(khat, 3)) * std::pow((1 + std::pow(khat, 2)), (-(7.0) / static_cast<double>(3.0)))) * directivity))) + 78.4);
        splhigh = (splhigh + (10.0 * std::log10((1 + (9.0 * std::pow(alphanoise, 2))))));
        sears = (1.0 / static_cast<double>((((twopi * kbar) / static_cast<double>(beta2)) + (1.0 / static_cast<double>((1.0 + ((2.4 * kbar) / static_cast<double>(beta2))))))));
        lfc = std::max<double>(aa_epsilon, ((((10 * sears) * mach) * std::pow(kbar, 2)) / static_cast<double>(beta2)));
        splti[(i - 1)] = (splhigh + (10.0 * log10aa((lfc / static_cast<double>((1 + lfc))))));
    }
    return splti;
}

Spectrum blunt(double alpstar, double c, double u, double theta, double phi, double l, double r, double h, double psi, const Parameters& p, double d99var2, double dstarvar1, double dstarvar2, double stallval) {
    double aterm = 0;
    double dbarh = 0;
    double deltap = 0;
    double dstarh = 0;
    double dstravg = 0;
    double dstrp = 0;
    double dstrs = 0;
    double eta = 0;
    double g4 = 0;
    double g50 = 0;
    double g514 = 0;
    double g5sum = 0;
    double hdstar = 0;
    double hdstarp = 0;
    int i = 0;
    double logval = 0;
    double m = 0;
    double rc = 0;
    double scale = 0;
    double stpeak = 0;
    double stppp = 0;
    Spectrum splblunt(p.freqlist.size(), 0.);
    Spectrum g5(p.freqlist.size(), 0.);
    m = (u / static_cast<double>(p.spdsound));
    rc = ((u * c) / static_cast<double>(p.kinvisc));
    if ((p.x_blmethod == x_blmethod_tables)) {
        deltap = d99var2;
        dstrs = dstarvar1;
        dstrp = dstarvar2;
    } else {
        std::tie(deltap, dstrs, dstrp) = thick(c, rc, alpstar, p, stallval);
    }
    dstravg = ((dstrs + dstrp) / static_cast<double>(2.0));
    hdstar = (h / static_cast<double>(dstravg));
    dstarh = (1.0 / static_cast<double>(hdstar));
    dbarh = directh_te(m, theta, phi);
    if ((dbarh <= 0)) {
        std::fill(splblunt.begin(), splblunt.end(), 0.0);
        return splblunt;
    }
    aterm = (0.212 - (0.0045 * psi));
    if ((hdstar >= 0.2)) {
        stpeak = (aterm / static_cast<double>(((1.0 + (0.235 * dstarh)) - (0.0132 * std::pow(dstarh, 2)))));
    } else {
        stpeak = (((0.1 * hdstar) + 0.095) - (0.00243 * psi));
    }
    if ((hdstar <= 5.0)) {
        g4 = (((17.5 * log10aa(hdstar)) + 157.5) - (1.114 * psi));
    } else {
        g4 = (169.7 - (1.114 * psi));
    }
    scale = (10.0 * log10aa(((((std::pow(m, 5.5) * h) * dbarh) * l) / static_cast<double>(std::pow(r, 2)))));
    g5sum = 0.0;
    for (i = 1; i < (p.freqlist.size() + 1); ++i) {
        stppp = ((p.freqlist[(i - 1)] * h) / static_cast<double>(u));
        eta = log10aa((stppp / static_cast<double>(stpeak)));
        g514 = g5comp(hdstar, eta);
        hdstarp = (((6.724 * std::pow(hdstar, 2)) - (4.019 * hdstar)) + 1.107);
        g50 = g5comp(hdstarp, eta);
        g5[(i - 1)] = (g50 + ((0.0714 * psi) * (g514 - g50)));
        if ((g5[(i - 1)] > 0.0)) {
            g5[(i - 1)] = 0.0;
        }
        g5sum = (std::pow(10, (g5[(i - 1)] / static_cast<double>(10))) + g5sum);
        if ((g5sum != 0)) {
            logval = std::max<double>(aa_epsilon, (1 / static_cast<double>(g5sum)));
        } else {
            logval = 1;
        }
        splblunt[(i - 1)] = (((g4 + g5[(i - 1)]) + scale) - (10 * std::log10(logval)));
    }
    return splblunt;
}

double g5comp(double hdstar, double eta) {
    double eta0 = 0;
    double g5 = 0;
    double k = 0;
    double m = 0;
    double mu = 0;
    if ((hdstar < 0.25)) {
        mu = 0.1211;
    } else {
        if ((hdstar <= 0.62)) {
            mu = ((-(0.2175) * hdstar) + 0.1755);
        } else {
            if ((hdstar < 1.15)) {
                mu = ((-(0.0308) * hdstar) + 0.0596);
            } else {
                mu = 0.0242;
            }
        }
    }
    if ((hdstar <= 0.02)) {
        m = 0.0;
    } else {
        if ((hdstar < 0.5)) {
            m = ((68.724 * hdstar) - 1.35);
        } else {
            if ((hdstar <= 0.62)) {
                m = ((308.475 * hdstar) - 121.23);
            } else {
                if ((hdstar <= 1.15)) {
                    m = ((224.811 * hdstar) - 69.354);
                } else {
                    if ((hdstar < 1.2)) {
                        m = ((1583.28 * hdstar) - 1631.592);
                    } else {
                        m = 268.344;
                    }
                }
            }
        }
    }
    m = std::max<double>(m, 0.0);
    eta0 = -(std::sqrt((((m * m) * std::pow(mu, 4)) / static_cast<double>((6.25 + (((m * m) * mu) * mu))))));
    if ((eta <= eta0)) {
        k = (((2.5 * std::sqrt((1.0 - std::pow((eta0 / static_cast<double>(mu)), 2)))) - 2.5) - (m * eta0));
        g5 = ((m * eta) + k);
    } else {
        if ((eta <= 0.0)) {
            g5 = ((2.5 * std::sqrt((1.0 - std::pow((eta / static_cast<double>(mu)), 2)))) - 2.5);
        } else {
            if ((eta <= 0.03615995)) {
                g5 = (std::sqrt((1.5625 - (1194.99 * std::pow(eta, 2)))) - 1.25);
            } else {
                g5 = ((-(155.543) * eta) + 4.375);
            }
        }
    }
    return g5;
}

double amin(double a) {
    double amina = 0;
    double x1 = 0;
    x1 = std::abs(a);
    if ((x1 <= 0.204)) {
        amina = (std::sqrt((67.552 - (886.788 * std::pow(x1, 2)))) - 8.219);
    } else {
        if ((x1 <= 0.244)) {
            amina = ((-(32.665) * x1) + 3.981);
        } else {
            amina = ((((-(142.795) * std::pow(x1, 3)) + (103.656 * std::pow(x1, 2))) - (57.757 * x1)) + 6.006);
        }
    }
    return amina;
}

double amax(double a) {
    double amaxa = 0;
    double x1 = 0;
    x1 = std::abs(a);
    if ((x1 <= 0.13)) {
        amaxa = (std::sqrt((67.552 - (886.788 * std::pow(x1, 2)))) - 8.219);
    } else {
        if ((x1 <= 0.321)) {
            amaxa = ((-(15.901) * x1) + 1.098);
        } else {
            amaxa = ((((-(4.669) * std::pow(x1, 3)) + (3.491 * std::pow(x1, 2))) - (16.699 * x1)) + 1.149);
        }
    }
    return amaxa;
}

double bmin(double b) {
    double bminb = 0;
    double x1 = 0;
    x1 = std::abs(b);
    if ((x1 <= 0.13)) {
        bminb = (std::sqrt((16.888 - (886.788 * std::pow(x1, 2)))) - 4.109);
    } else {
        if ((x1 <= 0.145)) {
            bminb = ((-(83.607) * x1) + 8.138);
        } else {
            bminb = ((((-(817.81) * std::pow(x1, 3)) + (355.21 * std::pow(x1, 2))) - (135.024 * x1)) + 10.619);
        }
    }
    return bminb;
}

double bmax(double b) {
    double bmaxb = 0;
    double x1 = 0;
    x1 = std::abs(b);
    if ((x1 <= 0.1)) {
        bmaxb = (std::sqrt((16.888 - (886.788 * std::pow(x1, 2)))) - 4.109);
    } else {
        if ((x1 <= 0.187)) {
            bmaxb = ((-(31.313) * x1) + 1.854);
        } else {
            bmaxb = ((((-(80.541) * std::pow(x1, 3)) + (44.174 * std::pow(x1, 2))) - (39.381 * x1)) + 2.344);
        }
    }
    return bmaxb;
}

double a0comp(double rc) {
    double a0 = 0;
    if ((rc < 95200.0)) {
        a0 = 0.57;
    } else {
        if ((rc < 857000.0)) {
            a0 = ((-(9.57e-13) * std::pow((rc - 857000.0), 2)) + 1.13);
        } else {
            a0 = 1.13;
        }
    }
    return a0;
}

std::tuple<double, double, double> thick(double c, double rc, double alpstar, const Parameters& p, double stallval) {
    double delta0 = 0;
    double deltap = 0;
    double dstr0 = 0;
    double dstrp = 0;
    double dstrs = 0;
    double logrc = 0;
    logrc = log10aa(rc);
    delta0 = (std::pow(10.0, ((1.6569 - (0.9045 * logrc)) + (0.0596 * std::pow(logrc, 2)))) * c);
    if ((p.itrip != itrip_none)) {
        delta0 = (std::pow(10.0, ((1.892 - (0.9045 * logrc)) + (0.0596 * std::pow(logrc, 2)))) * c);
    }
    if ((p.itrip == itrip_light)) {
        delta0 = (0.6 * delta0);
    }
    deltap = (std::pow(10.0, ((-(0.04175) * alpstar) + (0.00106 * std::pow(alpstar, 2)))) * delta0);
    if ((p.itrip != itrip_none)) {
        if ((rc <= 300000.0)) {
            dstr0 = ((0.0601 * std::pow(rc, -(0.114))) * c);
        } else {
            dstr0 = (std::pow(10.0, ((3.411 - (1.5397 * logrc)) + (0.1059 * std::pow(logrc, 2)))) * c);
        }
        if ((p.itrip == itrip_light)) {
            dstr0 = (dstr0 * 0.6);
        }
    } else {
        dstr0 = (std::pow(10.0, ((3.0187 - (1.5397 * logrc)) + (0.1059 * std::pow(logrc, 2)))) * c);
    }
    dstrp = (std::pow(10.0, ((-(0.0432) * alpstar) + (0.00113 * std::pow(alpstar, 2)))) * dstr0);
    if ((p.itrip == itrip_heavy)) {
        if ((alpstar <= 5.0)) {
            dstrs = (std::pow(10.0, (0.0679 * alpstar)) * dstr0);
        } else {
            if ((alpstar <= stallval)) {
                dstrs = ((0.381 * std::pow(10.0, (0.1516 * alpstar))) * dstr0);
            } else {
                dstrs = ((14.296 * std::pow(10.0, (0.0258 * alpstar))) * dstr0);
            }
        }
    } else {
        if ((alpstar <= 7.5)) {
            dstrs = (std::pow(10.0, (0.0679 * alpstar)) * dstr0);
        } else {
            if ((alpstar <= stallval)) {
                dstrs = ((0.0162 * std::pow(10.0, (0.3066 * alpstar))) * dstr0);
            } else {
                dstrs = ((52.42 * std::pow(10.0, (0.0258 * alpstar))) * dstr0);
            }
        }
    }
    return {deltap, dstrs, dstrp};
}

double directh_te(double m, double theta, double phi) {
    double dbar = 0;
    double degrad = 0;
    double mc = 0;
    double phir = 0;
    double thetar = 0;
    degrad = 0.017453;
    mc = (0.8 * m);
    thetar = (theta * degrad);
    phir = (phi * degrad);
    dbar = (((2.0 * std::pow(std::sin((thetar / static_cast<double>(2.0))), 2)) * std::pow(std::sin(phir), 2)) / static_cast<double>(((1.0 + (m * std::cos(thetar))) * std::pow((1.0 + ((m - mc) * std::cos(thetar))), 2))));
    return dbar;
}

double directh_le(double m, double theta, double phi) {
    double dbar = 0;
    double degrad = 0;
    double phir = 0;
    double thetar = 0;
    degrad = 0.017453;
    thetar = (theta * degrad);
    phir = (phi * degrad);
    dbar = (((2.0 * std::pow(std::cos((thetar / static_cast<double>(2.0))), 2)) * std::pow(std::sin(phir), 2)) / static_cast<double>(std::pow((1.0 + (m * std::cos(thetar))), 3)));
    return dbar;
}

double directl(double m, double theta, double phi) {
    double dbar = 0;
    double degrad = 0;
    double mc = 0;
    double phir = 0;
    double thetar = 0;
    degrad = 0.017453;
    mc = (0.8 * m);
    thetar = (theta * degrad);
    phir = (phi * degrad);
    dbar = (std::pow((std::sin(thetar) * std::sin(phir)), 2) / static_cast<double>(std::pow((1.0 + (m * std::cos(thetar))), 4)));
    return dbar;
}

Spectrum simple_guidati(double u, double chord, double thick_10p, double thick_1p, const Parameters& p) {
    double const1 = 0;
    double const2 = 0;
    int loop1 = 0;
    double slope = 0;
    double ti_param = 0;
    Spectrum splti(p.freqlist.size(), 0.);
    ti_param = (thick_1p + thick_10p);
    slope = ((1.123 * ti_param) + ((5.317 * ti_param) * ti_param));
    const1 = (((-(slope) * twopi) * chord) / static_cast<double>(u));
    const2 = (-(slope) * 5.0);
    for (loop1 = 1; loop1 < (p.freqlist.size() + 1); ++loop1) {
        splti[(loop1 - 1)] = ((const1 * p.freqlist[(loop1 - 1)]) + const2);
    }
    return splti;
}
}
