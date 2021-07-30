from random_variate_generation import rvg
import matplotlib .pyplot as plt

def unif():
    rv_maker = rvg(seed=101892)
    for i in range (1000000):
        rv_maker.unif(min=0, max=1, store=True)

    plt.hist(rv_maker.rn_storage["uniform"], rwidth=0.95)
    plt.show()

def normal():
    rv_maker = rvg(seed=101892)
    for i in range(1000000):
        rv_maker.norm(store=True)

    plt.hist(rv_maker.rn_storage["normal"], bins=40, rwidth=0.95)
    plt.show()

def bernoulli():
    rv_maker = rvg(seed=101892)
    for i in range(1000000):
        rv_maker.bern(p=.7, store=True)

    plt.hist(rv_maker.rn_storage["bernoulli"], bins=2, rwidth=0.95)
    plt.show()

def geometric():
    rv_maker = rvg(seed=101892)
    for i in range(1000000):
        rv_maker.geom(p=.5, store=True)

    plt.hist(rv_maker.rn_storage["geometric"], range=(0,12), rwidth=0.95)
    plt.show()

def exponential():
    rv_maker = rvg(seed=101892)
    for i in range(1000000):
        rv_maker.exp(lm=2, store=True)

    plt.hist(rv_maker.rn_storage["exponential"], rwidth=0.95)
    plt.show()

def gamma():
    rv_maker = rvg(seed=101892)
    for i in range(1000000):
        rv_maker.gamma(alpha=2, beta=2, store=True)

    plt.hist(rv_maker.rn_storage["gamma"], rwidth=0.95, bins=20)
    plt.show()


def weibull():
    rv_maker = rvg(seed=101892)
    for i in range(1000000):
        rv_maker.weibull(alpha=2, beta=2, store=True)

    plt.hist(rv_maker.rn_storage["weibull"], rwidth=0.95, bins=20)
    plt.show()

def poisson():
    rv_maker = rvg(seed=101892)
    for i in range(1000000):
        rv_maker.pois(lm=2, store=True)

    plt.hist(rv_maker.rn_storage["poisson"],range=(0,10), bins=10, rwidth=0.95)
    plt.show()

# def cauchy():
#     rv_maker = rvg(seed=101892)
#     for i in range(1000000):
#         rv_maker.cauchy(store=True)
#
#     plt.hist(rv_maker.rn_storage["cauchy"], rwidth=0.95)
#     plt.show()
#
# Did not have time to perfect this one :/

def triangle():
    rv_maker = rvg(seed=101892)
    for i in range(1000000):
        rv_maker.tri(mint=-10, maxt=10, store=True)

    plt.hist(rv_maker.rn_storage["triangular"], rwidth=0.95)
    plt.show()

def erlang():
    rv_maker = rvg(seed=101892)
    for i in range(1000000):
        rv_maker.erlang(lm=4, n=3, store=True)

    plt.hist(rv_maker.rn_storage["erlang"], rwidth=0.95)
    plt.show()

def gof(n):
    rv_maker = rvg(seed=101892)
    ao, ae = rv_maker.chi_squared_gof(n)

# unif()
# normal()
# bernoulli()
# geometric()
# exponential()
# gamma()
# weibull()
poisson()
# cauchy()
# triangle()
# erlang()

# gof(10000)