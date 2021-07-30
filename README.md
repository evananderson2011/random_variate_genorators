# random_variate_genorators

random_variate_generation.py is the importable class, used to create an object for RVG
rvg_test.py are sample tests that I have used to confirm the distributions being created are as expected, and that the chi squared GOF test is passed.

Random Variate Generation
The importable class "rvg" contained within this file is initialized with a seed value (for reproducability). After creation of the object (samples show in rvg_test), various random number based on given probabilities can be created. The core of this class is the random number generator "self.di()" that uses the desert island linear congruent we discussed in class. To ensure the generated values are uniform, the routine self.gof(n) can be run on any number of samples, splitting the results into 5 bins and conducting a 99% confident chi squared GOF test.

Using these randomly generated numbers, the class can generate random numbers sampled form a wide variety of distributions. For a given class, any time a random number is generated, the x_i value is incremented. This is important to keep in mind for reproducability. The x values are stored in self.x_list so can be values can be referenced more easily in the future. All variate generators for a given distribution have the option to be stored in the self.rn_storage disctionary. All routines can store these values by setting the optional parameter "store" to store=True.

Desert Island Generator - self.di():
This is the core RNG based on the LCG described in class. This LCG is the core of the uniform generator, and thus the backbone of random variate generation for all other distributions. It does not require inputs, and increments the x_i value each time a random variate is generated

Uniform:
A uniform distribution with any given minimum and maximum value. The values default to the standard uniform(0,1), but any minimum and maximum can be utilized. It returns a single value.

Normal:
A normal distribution random number generator, this routine utilizes the Box-Muller method to produce a value from a normal distrubtion based on the input mew and sigma values. It uses the formula X = mew + SQRT(sigma)*Z to transform the standard normal distribution norm(0,1) into the desired distribution before returning a single value. This routine defaults to producing the standard normal(0,1) distribution if no inputs are provided.

bernoulli:
A number randomly generated based on the bernoulli distribution, returning a single value (0 or 1) based on the required input p value. This routine returns 0 if the uniform generated is < 1 - p, otherwise returns 0

geometric:
A number randomly generated based on the inverse transform of the geometric distribution described in class. The required p value input is used to determine the output value.

exponential:
A randomly generated number based on the exponetial distrbution, derived from a uniform using the inverse transform described in class, with the required input lambda. This routine returns a single value.

gamma:
A randomly generated number based on the gamma distrbution, derived from a uniform using the inverse transform described in class, with the required inputs alpha and beta to describe the shape and scale respectively. This routine returns a single value.

weibull:
A randomly generated number based on the weibull distrbution, derived from a uniform using the inverse transform described in class, with the required inputs alpha and beta to describe the shape and scale respectively. This routine returns a single value.

poisson:
A randomly generated number based on the poisson distrbution, derived using mutliple uniforms and a variation of the acceptance-rejection method. The lambda is a required input, and this routine returns a single value.

triangle:
A randomly generated number based on the triangular distribution. The required input maxt and mint are the tri(max, min) values respectively to generate the distribution. This distributions assumes a center value as (mint + maxt)/2. The values are derived by generating two uniforms and adding them (the compostion method) to generate the number. This routine returns a single value.

erlang:
A randomly generated number based on the erlang_n distrbution, derived from a uniform using the inverse transform described in class, with the required inputs lambda and n. This routine returns a single value.


For user confidence, a routine to test the uniformity of the random numbers, self.chi_squared_gof(n) is included in the class as well. This routine samples n uniform values and usesthe chi squared goodness of fit test split into 5 bins to test the over all uniformity of the randomly generated uniforms.
