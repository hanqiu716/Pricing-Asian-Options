# Pricing-Asian-Options
Apply MC and QMC to price Asian options


Monte Carlo simulation (MC) is an approach that is widely used in high-dimensional
numerical integration and one of its main financial applications is option pricing. The
aim of this thesis is to evaluate the price of an Asian option using standard Monte
Carlo method and Quasi-Monte Carlo method (QMC) respectively. Since QMC's
convergence rate is determined by nominal problem dimension, the convergence rate
of QMC increases as the problem dimension increases, which limits the performance
of QMC in high dimensions. Hence, in this thesis, we also consider several techniques
which are proposed to capture the effective dimensions and improve the efficiency of
QMC in high-dimensional situations. The techniques include principal component
analysis (PCA) and Kronecker product approximation (KPA) and they are applied
for both constant and time-dependent volatilities. Finally, we conduct numerical
experiments and compare the precision and computational time between Quasi-Monte
Carlo and Monte Carlo methods.
