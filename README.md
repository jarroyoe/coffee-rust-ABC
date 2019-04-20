# coffee-rust-ABC
Approximate Bayesian Computation algorithm designed to run the Parameter estimation and uncertainty effects section in the research paper:
https://www.sciencedirect.com/science/article/pii/S0025556417306867

This project involved spatial simulations using NetLogo, whose spatial dynamics system was considered better for this research compared to that of R and the communication between R and NetLogo is relatively simple to implement.

In order to perform parameter estimation, we employed the EasyABC R Package applying the Beaumont algorithm depicted in https://academic.oup.com/biomet/article-abstract/96/4/983/220502

However, the connection and communication between R and NetLogo is computationally slow. To compensate this, we developed a parallel computation system which generates independent connections between R and NetLogo on each core of the processor.
This allowed the system to reduce its running time from the order of weeks to the order of days in a computer with 16GB of RAM and an Intel Core i7-8550U processor.
