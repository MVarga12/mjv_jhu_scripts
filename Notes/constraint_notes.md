# Gromacs Constraints

## Restraints to fix angles

  - Not harmonic
  - Can use one or two pairs of atoms; one pair will fix pitch or yaw, two pairs will fix both
  - Form is similar to fixing a dihedral:
    - one pair:
      $$
          V_{\text{ar}}(\mathbf{r}_i, \mathbf{r}_j) = k_{ar}(1-\cos(n(\theta - \theta_0)))
      $$

      where

      $$
          \theta = \arccos \left( \frac{\mathbf{r}_j - \mathbf{r}_i}{||\mathbf{r}_j - \mathbf{r}_i||} \cdot \left( \begin{smallmatrix} 0 \\ 0 \\ 1 \end{smallmatrix} \right)\right)
      $$

    - two pair:
      $$
          V_{\text{ar}}(\mathbf{r}_i, \mathbf{r}_j, \mathbf{r}_k, \mathbf{r}_l) = k_{ar}(1-\cos(n(\theta - \theta_0)))
      $$

      where

      $$
          \theta = \arccos \left( \frac{\mathbf{r}_j - \mathbf{r}_i}{||\mathbf{r}_j - \mathbf{r}_i||} \cdot \frac{\mathbf{r}_j - \mathbf{r}_k}{||\mathbf{r}_j - \mathbf{r}_k||}\right)
      $$

    - $n$ is the multiplicity (2 doesn't distinguish between parallel and anti-parallel vectors). $\theta$ should be between 0 and 180 degrees for $n=1$ and between 0 and 90 degrees for $n=2$

# Matlab scripts
    ```matlab
    % INPUTS: Reads in the coordinates for one COM (com1: x, y, z), which is therefore 
    % a vector with three elements: xcom, ycom, zcom. and the other COM. as
    % well as the boxlength in each dimension
    % Calculates the PBC corrected distance between the COM, and then the angle
    % between the vector and the z-axis.
    % OUPUTS: PBC corrected vector between the COMs and angle with z-axis.
    % Vector will be based on com2-com1!!!! Otherwise the angle will be
    % pi-angle

    function[vec, theta]=vector_and_angle(com1, com2, xbox, ybox, zbox)
    
    dx=com2(1)-com1(1);
    dy=com2(2)-com1(2);
    dz=com2(3)-com1(3);
    
    %Correct for PBC
    dx=dx-xbox/2*round(dx/(xbox/2));
    dy=dy-ybox/2*round(dy/(ybox/2));
    dz=dz-zbox/2*round(dz/(zbox/2));
    vec=[dx, dy, dz];
    
    veclen=sqrt(vec*vec');
    theta_rad=acos(dz/veclen);
    theta=theta_rad*180/pi
    ```

# Gromacs command notes
  - `editconf -c -box` takes a `.gro` file and orients it so that the center of mass lies at the center of some xyz box.
    - example:
      ```bash
      gmx editconf -f some.gro -o some_oriented.gro -c -box 6.4 6.5 6.6
      ```
      orients `some.gro` such that its center of mass, denoted by `c`, is the center of the (PBC cell) box with dimensions 6.4 6.5 6.6
  - solvating systems involving membranes is difficult, since it tends to fill empty space between the lipids with water molecules. The author of the tutorials suggests making a local copy of `vdwradii.dat` and changing the van der Waals radius of carbon from 0.17 to 0.375 for use in solvation ONLY

