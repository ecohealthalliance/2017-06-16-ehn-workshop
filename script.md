# Simulating Epidemiological Models with R
Noam Ross  
2017-06-05  

_This material is adapted from
[lessons from the 2012 EEID workshop](https://ms.mcmaster.ca/~bolker/eeid/)
by Helen J. Wearing, John M. Drake & Aaron A. King_

# Introduction - logistics

-  I'm Noam. EHA diease ecologist, background in theoretical ecology. Do a mix of 
   applied and theoretical modeling, focus on dynamics. 
-  Here today to talk about epidemiological modeling.
-  It's a big topic, and it's Friday at the end of a long week.
-  So I'll aim to riff on some theoretical concepts.  Then we'll try to dive
   in as fast as possible into hands-on work.
-  For our hands-on work we are going to use the R programming language.
   It's not he easiest thing to use but it gives you the tools to adapt what
   you learn to your own work.
-  Now, many of you took my poll on what you knew about certain modeling
   and programming topics. About half of you had some familiarity with R.
   If you have, please raise your hand.
-  I'd like you to pair up with someone who does not have their hand up.  We're
   practice "pair programming". Much of the lesson will be me live-coding on
   the screen and you following along.  Only one of you should follow-it doesn't
   matter which.
-  Secondly, you each have a green and red post-it. If you are having trouble,
   if you've fallen behind, put up a red stickie like this.
-  Laurie, who has some experience with this work, has agreed to act as my TA,
   will dive in to help.  When we do exercises, we'll use red to ask for help
   and green to indicate you're done with a problem.
-  Also, everything I type up here will be available to you as we go. If
   you go to bit.ly/ehn2017, you can refresh the browser to update everything
   that I'm doing.
   
# So, let's dive in!

In this lesson, I introduce some model structures that are useful for
representing the **temporal dynamics** of disease transmission. These are simple
models for representing how change in time represents interactions.   
the models we look at here are fundamental and relatively general and
therefore readily extended for your own purposes in the future. Third, we
introduce some numerical tools that are useful for studying epidemiological
systems. 

# Why do these types of models?

-   These models are abstract representations of simplified systems.  As the
    saying goes, "All models are wrong, some models are useful."  So, how
    are these very simple models _useful_?

-  Models of this sort tend not to be very _predictive_.   For models, you
   have to balance generality (representing transferable mechanisms),
   accuracy (detailed representation of the system), and preditiveness. 
   
-   These simple models are _general_, and this is their advantage.

-   They are readily extended and modified for your own purposes in the future
-   Their simplicity lets us explore the implications of a small number of
        key mechanisms.  The insights we gain from studying a model help us understand
        of a set of interactions gives rise to certain behaviors, and then we
        us these to make inferences about how real-world systems proceed. 
-  We can apply our qualitative insights from them across multple systems.
-  They are also good building blocks to extend to develop models for your
own purposes in the future.
        

# Chain binomial model

We begin by developing an intuitive understanding of the mechanics of the
transmission process by considering a simple model of an epidemic:
the **chain binomial** model.

This is one of the very first epidemic models. It is most commonly attributed
Lowell Reed and Wade Hamton Frost's work in the 1920s, and hence called the Reed-Frost model, but
a formulation of it was actually published in the 1880s by a physician-mathematician
names Piotr Dmitirievich En'ko, in Russia.

What are some properties of this model?

- _Compartmental_ model
- _Discrete_ in time and in population
- _Non-overlapping_ generations
- _Stochastic_

We designate the number of susceptible individuals at a given time as $I_t$, 
and the infected individuals as $S_t$

This is a very simple model with types of organisms - susceptible and infected,
and two parameters, the contact rate $\beta$, and the number of susceptibles,
$S_0$.
df
This model stipulates that the epidemic evolves according to discrete
generations. In each generation, new infections are binomially distributed
with the number of trials equal to the number of susceptibles, $S_t$, and
probability of infection, $p=1-\exp(-\beta - I_t)$. In probability notation:

$$I_{t+1} \sim \mathrm{binom}(S_{t},1-\exp(-\beta\,I_{t}))$$

Susceptibles are then depleted by the number of these infections

$$S_{t+1} = S_{t}-I_{t+1}$$

Some of you may have taken statistics or probability courses and learned that the binomial random variable is
the number of independent "successes" in a sequence of weighted coin
tosses. The analogy here is that we toss a weighted coin (with probability of
heads $p=1-\exp(-\beta\,I_t)$) for each susceptible individual in the
population. If the weighted coin does in fact come up heads then the
susceptible individual becomes infected. Otherwise, it stays susceptible and
we move on to the next susceptible individual.


```r
state_0 <- c(S=2000, I=1) # do three things here.  I create an object (state_0), I make an object with two variables in it, and
state_0
```

```
##    S    I 
## 2000    1
```

```r
my_parms <- c(beta=0.001)
my_parms 
```

```
##  beta 
## 0.001
```

Now we have to express that equation in code.  To do so we write a 
_function_.  A function is a reusable bit of code that takes something
as an input, and returns something as an output.  We'll write a function
that executes one step of this binomial chain.


```r
binomial_chain_step = function(state, parameters) {
  S = state["S"]
  I = state["I"]
  beta = parameters["beta"]
  I_next = rbinom(n = 1, size = S, prob = 1 - exp(-beta*I))
  S_next = S - I_next
  state_next = c(S = S_next, I = I_next)
  return(state_next)
}
```

Now let's run this function, using our state at parameters as inputs.


```r
binomial_chain_step(state_0, my_parms)
```

```
##  S.S    I 
## 1999    1
```

We can do this several times and we see the randomness inherent in this model.

Now we write another function that runs this many times a row in a _loop_,
and store the results each time. It needs another input - how many steps to run?


```r
binomial_chain_simulate = function(state, parameters, n_steps) {
  output = matrix(nrow = n_steps + 1, ncol = 3)
  output[1,1] = 0
  output[1,2] = state["S"]
  output[1,3] = state["I"]
  colnames(output) <- c("time", "S", "I")
  for (step in 1:n_steps) {
    output[step + 1, 1] = step
    output[step + 1, 2:3] = binomial_chain_step(output[step, 2:3], parameters)
  }
  return(output)
}
```


```r
binomial_chain_simulate(state_0, my_parms, 10)
```

```
##       time    S   I
##  [1,]    0 2000   1
##  [2,]    1 1999   1
##  [3,]    2 1998   1
##  [4,]    3 1996   2
##  [5,]    4 1992   4
##  [6,]    5 1983   9
##  [7,]    6 1967  16
##  [8,]    7 1937  30
##  [9,]    8 1881  56
## [10,]    9 1763 118
## [11,]   10 1568 195
```

```r
results <- binomial_chain_simulate(state_0, my_parms, 100)
```

We can plot these results. 



```r
plot(I~time, data=results, type="l")
```

![](script_files/figure-html/unnamed-chunk-6-1.png)<!-- -->


We'll now specify some parameters, simulate the model a few times, and plot
the results.



```r
n_sims <- 100
n_step <- 30
results_all <- list()
for (k in 1:n_sims) {
  results_all[[k]] <- binomial_chain_simulate(state_0,my_parms, n_step)
}
```


```r
plot(c(0,20),c(0,400),type="n",xlab="time",ylab="I")
for (k in 1:n_sims) {
  lines(I~time, data = results_all[[k]], type="l")
}
```

![](script_files/figure-html/unnamed-chunk-8-1.png)<!-- -->

We can save this plot using the "Export" button on top of the Plot in RStudio.


>## Exercise
>
>Explore the dynamics of the system for different values of
>$\beta$, as well as different initial values of $S$ and $I$. 

Although simple, the chain binomial model captures some key properties of the real biological process:

-   demographic stochasticity - a type of process noise (to be contrasted with **measurement error**)
-   item categorical class variables ($S$ and $I$ are integer-valued)

# Two-species model

We now extend this model to a *two species* version to look at spillover scenarios.
Now we have to classes each of susceptible and infected individuals: $S1, I1, S2, I2$


```r
state2_0 = c(S1 = 2000, I1 = 1, S2 = 2000, I2 = 0)
parameters2 = c(beta1 = 0.001, beta2 = 0.001, beta12 = 0.0001)
```




```r
binomial_chain_step2 = function(state, parameters) {
  S1 = state["S1"]
  S2 = state["S2"]
  I1 = state["I1"]
  I2 = state["I2"]
  beta1 = parameters["beta1"]
  beta2 = parameters["beta2"]
  beta12 = parameters["beta12"]
  I1_next = rbinom(n = 1, size = S1, prob= 1 -exp(-beta1 * I1))
  S1_next = S1 - I1_next
  I2_next = rbinom(n = 1, size = S2, prob = 1 - exp(-(beta2 * I2 + beta12 * I1)))
  S2_next = S2 - I2_next
  return(c(S1_next, I1_next, S2_next, I2_next))
}
```


```r
binomial_chain_step2(state2_0, parameters2)
```

```
##   S1        S2      
## 1999    1 1999    1
```


Then we put these together in sequence to simulate an entire epidemic:



```r
binomial_chain_simulate2 = function(state, parameters, n_steps) {
  output = matrix(nrow = n_steps + 1, ncol = 5)
  colnames(output) <- c("time", "S1", "I1", "S2", "I2")
  output[1,1] <- 0
  output[1, 2:5] <- state
  for (step in 1:n_steps) {
    output[step + 1, 1] <- step
    output[step +1, 2:5] <- binomial_chain_step2(output[step, 2:5], parameters)
  }
  return(output)
}
```


```r
binomial_chain_simulate2(state2_0, parameters2, 100)
```

```
##        time   S1  I1   S2  I2
##   [1,]    0 2000   1 2000   0
##   [2,]    1 1999   1 1999   1
##   [3,]    2 1996   3 1998   1
##   [4,]    3 1994   2 1995   3
##   [5,]    4 1990   4 1989   6
##   [6,]    5 1983   7 1966  23
##   [7,]    6 1966  17 1923  43
##   [8,]    7 1938  28 1823 100
##   [9,]    8 1879  59 1632 191
##  [10,]    9 1787  92 1327 305
##  [11,]   10 1629 158  959 368
##  [12,]   11 1406 223  659 300
##  [13,]   12 1124 282  456 203
##  [14,]   13  857 267  368  88
##  [15,]   14  659 198  321  47
##  [16,]   15  543 116  295  26
##  [17,]   16  486  57  282  13
##  [18,]   17  459  27  279   3
##  [19,]   18  449  10  278   1
##  [20,]   19  444   5  277   1
##  [21,]   20  442   2  277   0
##  [22,]   21  442   0  277   0
##  [23,]   22  442   0  277   0
##  [24,]   23  442   0  277   0
##  [25,]   24  442   0  277   0
##  [26,]   25  442   0  277   0
##  [27,]   26  442   0  277   0
##  [28,]   27  442   0  277   0
##  [29,]   28  442   0  277   0
##  [30,]   29  442   0  277   0
##  [31,]   30  442   0  277   0
##  [32,]   31  442   0  277   0
##  [33,]   32  442   0  277   0
##  [34,]   33  442   0  277   0
##  [35,]   34  442   0  277   0
##  [36,]   35  442   0  277   0
##  [37,]   36  442   0  277   0
##  [38,]   37  442   0  277   0
##  [39,]   38  442   0  277   0
##  [40,]   39  442   0  277   0
##  [41,]   40  442   0  277   0
##  [42,]   41  442   0  277   0
##  [43,]   42  442   0  277   0
##  [44,]   43  442   0  277   0
##  [45,]   44  442   0  277   0
##  [46,]   45  442   0  277   0
##  [47,]   46  442   0  277   0
##  [48,]   47  442   0  277   0
##  [49,]   48  442   0  277   0
##  [50,]   49  442   0  277   0
##  [51,]   50  442   0  277   0
##  [52,]   51  442   0  277   0
##  [53,]   52  442   0  277   0
##  [54,]   53  442   0  277   0
##  [55,]   54  442   0  277   0
##  [56,]   55  442   0  277   0
##  [57,]   56  442   0  277   0
##  [58,]   57  442   0  277   0
##  [59,]   58  442   0  277   0
##  [60,]   59  442   0  277   0
##  [61,]   60  442   0  277   0
##  [62,]   61  442   0  277   0
##  [63,]   62  442   0  277   0
##  [64,]   63  442   0  277   0
##  [65,]   64  442   0  277   0
##  [66,]   65  442   0  277   0
##  [67,]   66  442   0  277   0
##  [68,]   67  442   0  277   0
##  [69,]   68  442   0  277   0
##  [70,]   69  442   0  277   0
##  [71,]   70  442   0  277   0
##  [72,]   71  442   0  277   0
##  [73,]   72  442   0  277   0
##  [74,]   73  442   0  277   0
##  [75,]   74  442   0  277   0
##  [76,]   75  442   0  277   0
##  [77,]   76  442   0  277   0
##  [78,]   77  442   0  277   0
##  [79,]   78  442   0  277   0
##  [80,]   79  442   0  277   0
##  [81,]   80  442   0  277   0
##  [82,]   81  442   0  277   0
##  [83,]   82  442   0  277   0
##  [84,]   83  442   0  277   0
##  [85,]   84  442   0  277   0
##  [86,]   85  442   0  277   0
##  [87,]   86  442   0  277   0
##  [88,]   87  442   0  277   0
##  [89,]   88  442   0  277   0
##  [90,]   89  442   0  277   0
##  [91,]   90  442   0  277   0
##  [92,]   91  442   0  277   0
##  [93,]   92  442   0  277   0
##  [94,]   93  442   0  277   0
##  [95,]   94  442   0  277   0
##  [96,]   95  442   0  277   0
##  [97,]   96  442   0  277   0
##  [98,]   97  442   0  277   0
##  [99,]   98  442   0  277   0
## [100,]   99  442   0  277   0
## [101,]  100  442   0  277   0
```

```r
results2 <- binomial_chain_simulate2(state2_0, parameters2, 100)
```


```r
plot(I1~time, data = results2, type="l", col="blue")
lines(I2~time, data= results2, type="l", col="red")
```

![](script_files/figure-html/unnamed-chunk-14-1.png)<!-- -->



```r
n_sims = 1000
n_steps = 30

state2_0 = c(S1 = 2000, I1 = 1, S2 = 2000, I2 = 0)
parameters2 = c(beta1 = 0.0005, beta2 = 0.005, beta12 = 0.00001)
results_all2 = list()
for (k in 1:n_sims) {
  results_all2[[k]] = binomial_chain_simulate2(state2_0, parameters2, n_steps)
}
```


```r
plot(c(0,30), c(0,1500), type="n", xlab="time", ylab="I")
for (k in 1:n_sims) {
  lines(I1~time, data=results_all2[[k]], type="l", col="blue")
  lines(I2~time, data=results_all2[[k]], type="l", col="red")
}
```

![](script_files/figure-html/unnamed-chunk-16-1.png)<!-- -->

>## Challenge
>
> 1. Simulate a disease, such as rabies, where spillover is common from the reservoir host to 
>    the spillover host, but there is little spread in the spillover hosts.
> 2.  Simulate a disease, such as influenza, where spillover from wild hosts
>     is extremely rare, but the disease is highly contagious in spillover
>     hosts
> 3. Make the disease in (2) less contagious in the wild host.  What type of
>    patterns do we see?

----
*This is as far as we at the EcoHealth Net Workshop*

# SIR Models

However, the chain binomial, like all models, is an approximation. One large assumption that it makes is that the generations are perfectly synchronized. For some diseases, this may not be such a bad approximation; for others, it might very well be.

There is an implicit parameter here - that the period of time that the diease
lasts is equal to 1.  Both your model structure and parameters can contain assumptions.

Let's have a look at what can be done with models that don't make this assumption, i.e. generations of infection are not synchronized. In fact, for now, we'll take it one step further and assume that the change in the number of susceptible and infectious individuals in the population happens continuo=ntain assumptions.

The simplest place to start is with the classical $SIR$ model. This model divides the host population into three classes with respect to their infection status:
individuals are either Susceptible, Infected (and Infectious), or Recovered.



The model simply keeps track of how many individuals are in each class: individuals that leave one class must enter another.
The only exceptions, of course, are births and deaths. 

The SIR model makes some assumptions: we're going to assume we have a large (technically infinitely large) population in which the effects of demographic stochasticity become negligible.  Again,
this can be a useful simplifying assumption.  There are options for stochastic,
overlapping generations models, of course, too.

The state variables change according to a
system of differential equations:

$$\begin{aligned}
\frac{dS}{dt} &= B-\lambda(I,t)\,S-\mu\,S\\
\frac{dI}{dt} &= \lambda(I,t)\,S-\gamma\,I-\mu\,I\\
\frac{dR}{dt} &= \gamma\,I-\mu\,R\\
\end{aligned}$$

Here, $B$ is the crude birth rate, $\mu$ is the per capita death rate, $N$ is the host population size, and $\gamma$ the recovery
rate.  The term that makes this model interesting (and nonlinear) is the **force-of-infection**, represented by the function
$\lambda(I,t)$.  We'll assume that it has the so-called
**frequency-dependent** form

$$\lambda(I,t) = \beta(t)\,\frac{I}{N}$$

So that the risk of infection faced by a susceptible individual is proportional to the fraction of the population that is infectious.  Notice that we allow for the possibility of a contact rate, $\beta$, that varies in time. In this model, $S$, $I$, and $R$ may be interpreted either as proportions of the population (if $N=1$) or abundances (if $N>1$).

Like many epidemiological models, one can't solve the $SIR$ equations
explicitly.  Rather, to find the trajectory of a continuous-time model such as the $SIR$, we must integrate those ordinary differential equations (ODEs) numerically.  What we mean by this is that we use a computer algorithm to approximate the solution.  In general, this can be a tricky business. Fortunately, this is a well studied problem in numerical analysis and (when the equations are smooth, well-behaved functions of a relatively small number of variables) standard numerical integration schemes are available to approximate the integral with arbitrary precision.  Particularly, \R\ has very sophisticated ODE solving capabilities in the package \code{deSolve}. To use these algorithms we first load the package:



```r
# When running the script on your computer, you may need to to run
# install.packages("deSolve")
# to install this package before proceeding
library(deSolve)
```

The ODE solver needs to know the right-hand sides of the ODE.
We give it this information as a function:


```r
sir_diff_eqs <- function (time, state, parameters) {
  ## first extract the state variables
  S <- state["S"]
  I <- state["I"]
  R <- state["R"]
  N <- S+I+R
  ## now extract the parameters
  beta <-  parameters["beta"]
  gamma <- parameters["gamma"]
  mu <-    parameters["mu"]
  B <-     parameters["B"]
  ## now code the model equations
  dSdt <- B-beta*S*I/N-mu*S
  dIdt <- beta*S*I/N-(mu+gamma)*I
  dRdt <- gamma*I-mu*R
  ## combine results into a single vector
  derivs <- list(c(dSdt,dIdt,dRdt))
  ## return result as a list!
  return(derivs)   
}
```


Notice that in this case we've assumed $\beta$ is constant.



We'll also write a function to calculate $R_0$.



```r
R0 <- function(parameters) {
  R0 = parameters["beta"] /
         (parameters["mu"]+parameters["gamma"])
}
```


We'll now define the times at which we want solutions, assign some
values to the parameters, and specify the **initial conditions**,
**i.e.**, the values of the state variables $S$, $I$, and $R$ at the
beginning of the simulation:



```r
times <- seq(0,30,by=1/120)
parameters  <- c(B=1/70,mu=1/70,N=1,beta=400,gamma=365/14)
state_0 <- c(S=1-0.001-0.9,I=0.001,R=0.9)
```


Now we can simulate a model trajectory with the \code{ode} command:



```r
out <- ode(state_0,times,sir_diff_eqs,parameters)
```


and plot the results



```r
plot(R~time,data=out,type='l', ylim=c(0,1), ylab="Population")
lines(I~time,data=out,type='l',col='red'); par(new=TRUE)
lines(S~time,data=out,type='l',  col='blue')
```

![](script_files/figure-html/unnamed-chunk-22-1.png)<!-- -->

The "cycling" is even more apparent if we look at comparisons of 


```r
plot(I~S,data=out) #,type='b',log='xy',yaxt='n',xlab='S',cex=0.5)
```

![](script_files/figure-html/unnamed-chunk-23-1.png)<!-- -->

> ## Challenges
>
> Explore the dynamics of the system for different values of the
> $\beta$ and $B$ parameters by simulating and plotting trajectories
> as time series and in phase space (e.g., $I$ vs. $S$).  How the $\beta$, $B$,
> and $R0$ related to the type of trajectories you get?

> ## Big Challenges (Pick one, if there is time, do as group)
>  1.  What if you have a disease that has a latent period in the host before
>  it starts infecting other hosts?  Change the SIR model to an SEIR model
>  with four compartments (Susceptible, EXPOSED, Infectious, Recovered) and
>  plot the dynamics of all four.
>  2.  What if there are seasonal dynamics to the disease? Change the model
>  so that the value of either $\B$ cycles up and down, representing seasonality
>  of reproduction, or $\beta$ cycles up and down, representing seasonality in
>  contact.  (Hint: there are `sin()` and `cos()` functions, and `time` is 
>  a variable in your differential equation function already).  Plot several
>  parameterizations of this model and compare to the version already created.
