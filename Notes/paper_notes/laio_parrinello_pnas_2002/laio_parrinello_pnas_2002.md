# Notes on Laio & Parrinello, PNAS, 2002
## Summary

Paper which introduces metadynamics. MD and MC simulations, used in the context of situations in which the free energy
surface (FES) of the system has local minima separated by large free energy barriers, are limited by their short time
frame and sampling algorithms. As a result, simulations started in one minimum are not likely to move to the next
minimum, over the large free energy barrier.

## Methodology

- assume finite number of collective coordinates \(s_i,~i=1\), where \(n\) is a small number
  - we consider the dependence of the free energy, \(\mathcal{F}(s)\), on these collective coordinates
  - exploration of the FES is controlled by forces \(F^t_i = \delta \mathcal{F}/\delta s_i^t\)
  - note that \(t\) is an index of the states

- introduce ensemble of \(P\) replicas, each with constraints on the collective coordinates such that \(s_i = s_i^t\)
  and evolve separately at the same temperature
  - the \(P\) replicas are statistically independent --- can produce statistics faster in parallel than in serial
  - constraints imposed on the replicas by adding \(\sum_i = _{1,n} \lambda_i (s_i-s_i^t)\) to the Lagrangean
  - can obtain the derivative of the free energy, \(F_i^t = \langle \lambda_i \rangle\), relative to \(s_i^t\), as an
    average over time and the replicas

- use the forces to define coarse--grained dynamics 
  - defined from the discretized evolution equation:
  \[
  \sigma_i^{t+1} = \sigma_i^t + \delta\sigma \frac{\phi^t_i}{\phi^t}
  \]
  where \(\sigma_i^t = s_i^t/\Delta s_i\) are the scaled variables, \(\phi_i^t = F_i^t \Delta s_i\) are the scaled
  forces, \(\Delta s_i\) is the estimated size of the FES well in the direction \(s_i\), \(|\phi^t|\) is the modulus of
  the \(n\)-dimensional vector \(\phi_i^t\), and \(\delta\sigma\) is a dimensionless stepping parameter.
  - after the collective coordinates are updated via the above discretized evolution equation, a new ensemble of
    replicas with values \(\sigma_i^{t+1}\) is prepared, and new forces \(F_i^{t+1}\) are calculated for the next
    iteration
  - There is no dynamical continuity between the replicas, so one can use large values of \(\delta\sigma\)

- need to replace the forces with add history dependence

  \[ \sigma_i \rightarrow \sigma_i - \frac{\delta}{\delta\sigma_i} W \sum\limits{t\prime \leq t} \product\limits_i \]
