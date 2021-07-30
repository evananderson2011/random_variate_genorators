# This is a library for generating random variates built by @notevananderson

import numpy
import pandas
from scipy import stats as stats
from math import fmod as mod
from numpy import log as ln
import math


class rvg:
    def __init__(self, seed):

        if (seed < 1) or (seed > 2**31) or not isinstance(seed, type(int())):
            error = "The seed value must be between 1 and 2 to the 31st power and must be an integer"
            raise error

        self.x_0 = seed
        self.x_i = seed
        self.x_list = [].append(self.x_0)
        self.di_m = ((2 ** 31) - 1)

        self.rn_storage = {
            "bernoulli": [],
            "geometric": [],
            "exponential": [],
            "normal": [],
            "gamma": [],
            "weibull": [],
            # "cauchy": [],
            "poisson": [],
            "cauchy": [],
            "uniform": [],
            "triangular": [],
            "erlang": []
            }

    def di(self):
        # Desert island RN generator
        self.x_i = int(mod((16807 * self.x_i), self.di_m))
        self.x_list = [].append(self.x_i)

        return self.x_i

    def unif(self, min=0, max=1, store=False):
        # Using the desert island generator for all uniforms
        # rounding all uniforms to 12 places for storage and speed purposes
        u_i = self.di() / self.di_m
        u_i = (u_i * (max - min)) + min

        if store:
            self.rn_storage["uniform"].append(u_i)

        return u_i

    def norm(self, mew=0, sigma2=1, store=False):
        if sigma2 < 0:
            error = "Sigma squared cannot be less than zero"
            raise error

        # using box-muller method to generate:
        u_1 = self.unif()
        u_2 = self.unif()

        z_1 = math.sqrt(-2*ln(u_1)) * math.cos(2 * math.pi * u_2)  # This is normally distributed

        if mew != 0 or mew != 1:
            rv_i = mew + math.sqrt(sigma2) * z_1
        else:
            rv_i = z_1

        if store:
            self.rn_storage["normal"].append(rv_i)

        return rv_i

    def bern(self, p, store=False):
        if p < 0 or p > 1:
            error = "The p value for bernoulli distributions must be between zero and 1"
            raise error

        u_i = self.unif()

        if u_i < 1 - p:
            rv_i = 0
        else:
            rv_i = 1

        if store:
            self.rn_storage["bernoulli"].append(rv_i)

        return rv_i

    def geom(self, p, store=False):
        if p < 0 or p > 1:
            error = "The p value for geometric distributions must be between zero and 1"
            raise error

        rv_i = ln(self.unif())/ln(1-p)
        rv_i = math.ceil(rv_i)

        if store:
            self.rn_storage["geometric"].append(rv_i)

        return rv_i

    def exp(self, lm, store=False):
        rv_i = -(1/lm) * ln(self.unif())

        if store:
            self.rn_storage["exponential"].append(rv_i)

        return rv_i

    def gamma(self, alpha, beta, store=False):
        lm = alpha

        if beta < 1:
            b = (math.e + beta)/math.e
            while True:
                u_i = self.unif()
                W = b * u_i

                if W < 1:
                    Y = W ** (1/beta)
                    V = self.unif()
                    if V <= math.e ** (-1*Y):
                        rv_i = Y/lm
                        break

        if beta >= 1:
            a = ((2 * beta) - 1) ** (-1/2)
            b = beta - ln(4)
            c = beta + (a ** -1)
            d = 1 + ln(4.5)

            while True:
                u_1 = self.unif()
                u_2 = self.unif()

                V = a * ln(u_1/(1 - u_1))
                Y = beta * (math.e ** V)
                Z = (u_1 ** 2)*(u_2)
                W = b + (c*V)  - Y

                if W + d - (4.5*Z) >= 0:
                    rv_i = Y/lm
                    break
                elif W >= ln(Z):
                    rv_i = Y/lm
                    break

        if store:
            self.rn_storage["gamma"].append(rv_i)

        return rv_i

    def weibull(self, alpha, beta, store=False):
        rv_i = (1 / alpha) * ((-1 * ln(self.unif())) ** (1/beta))

        if store:
            self.rn_storage["weibull"].append(rv_i)

        return rv_i

    def pois(self, lm, store=False):
        # Functionally doing the discrete inverse transform of the c.d.f.
        x = 0
        bigF_x = 0

        while True:
            u_i = self.unif()
            f_x = ((math.e ** -lm) * (lm ** x))/(math.factorial(x))
            bigF_x += f_x
            if u_i < bigF_x:
                break
            x += 1

        rv_i = x
        if store:
            self.rn_storage["poisson"].append(rv_i)

        return rv_i

    # Cauchy distribution was not finished do to lack of time :/
    # def cauchy(self, store=False):
    #     rv_i = math.tan(2*math.pi*self.unif())
    #
    #     if store:
    #         self.rn_storage["cauchy"].append(rv_i)

    def tri(self, mint, maxt, store=False):

        u_1 = self.unif(min=mint/2, max=maxt/2)
        u_2 = self.unif(min=mint/2, max=maxt/2)

        rv_i = u_1 + u_2

        if store:
            self.rn_storage["triangular"].append(rv_i)

        return rv_i

    def erlang(self, lm, n=1, store=False):

        if n == 1:
            notice = "When n=1, erlang dist. is exponential dist, please use the self.exp() function"
            raise notice
        if n < 1:
            error = "n must be larger than 2. If n=1 is desired, use self.exp() instead()"
            raise error
        product = 1

        for i in range(n):
            product = product * self.unif()

        rv_i = (-1/lm)*ln(product)

        if store:
            self.rn_storage["erlang"].append(rv_i)

        return rv_i

    def chi_squared_gof(self, n):

        # Performing the xhi squared goodness of fit test for 1000 Unif(0.1) values, split in to 5 bins:
        # [0, .2), [.2, .4), [.4, .6), [.6, .8), [.8, 1]
        exp_vs_actual = {
            "exp":int(n/5),
            "one":0,
            "two": 0,
            "three": 0,
            "four": 0,
            "five": 0,
        }

        for i in range(n):
            val_i = self.unif()
            if val_i < .2:
                exp_vs_actual["one"] += 1
            elif val_i < .4:
                exp_vs_actual["two"] += 1
            elif val_i < .6:
                exp_vs_actual["three"] += 1
            elif val_i < .8:
                exp_vs_actual["four"] += 1
            else:
                exp_vs_actual["five"] += 1

        chi_2_obv = ((exp_vs_actual["one"] - exp_vs_actual["exp"])**2)/exp_vs_actual["exp"]
        chi_2_obv += ((exp_vs_actual["two"] - exp_vs_actual["exp"]) ** 2) / exp_vs_actual["exp"]
        chi_2_obv += ((exp_vs_actual["three"] - exp_vs_actual["exp"]) ** 2) / exp_vs_actual["exp"]
        chi_2_obv += ((exp_vs_actual["four"] - exp_vs_actual["exp"]) ** 2) / exp_vs_actual["exp"]
        chi_2_obv += ((exp_vs_actual["five"] - exp_vs_actual["exp"]) ** 2) / exp_vs_actual["exp"]

        chi_2_alpha = stats.chi2.ppf(.99, df=n-1)

        print('With {} observation, a chi squared value of {} was observed. \n'.format(n, round(chi_2_obv, 3)))
        if chi_2_obv < chi_2_alpha:
            print('This is less than the test statistic for alpha = 0.01 and k-1 observations, {}\nTherefore we fail to reject the null hypothesis, implying that the RNGs are uniform'.format(round(chi_2_alpha, 3)))
        elif chi_2_obv < chi_2_alpha:
            print('This is more than the test statistic for alpha = 0.05 and k-1 observations, {}\nTherefore we reject the null hypothesis that the RNGs are uniform, implying these randoms numbers are not uniform'.format(round(chi_2_alpha,3)))
        return chi_2_obv, chi_2_alpha



    # Helper to clear out storage if we want to use the same object
    def erase_storage(self):
        self.rn_storage = {
            "bernoulli": [],
            "geometric": [],
            "exponential": [],
            "normal": [],
            "gamma": [],
            "weibull": [],
            "conchy": [],
            "poisson": [],
            "cauchy": [],
            "uniform": [],
            "triangular": [],
            "erlang": []
        }







