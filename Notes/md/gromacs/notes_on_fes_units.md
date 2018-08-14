# Notes on Plumed FES

## Plumed Default Units
  - Units are agnostic of which MD package is being used
    - Energy = kJ/mol
    - Length = nm
    - Time = ps

### Units of Various Files
  - Angle constraint:
    - Using `vector_and_angle.py` and instructions laid out in `constraint_notes.pdf`, -ln(PDF)\*kT, with k = kJ/(mol\*K), units are in kJ/mol
  - Distance contstraints:
    - Using `calc_pos_restraint.py`, instructions laid out in `constraint_notes.pdf`, and k = kJ/(mol\*K), units are in kJ/mol, applied only to z-axis

## fes.dat example

## Dariush's conversion script

    % Creating plot for FES from PLUMED
    clear 

    data_filename = 'fes.dat';

    b_table = readtable(data_filename); b = table2array(b_table);

    my_legend1 = 'Distance 1';
    my_legend2 = 'Distance 2';

    my_title1 = 'Distance vs Time';
    my_title2 = 'Free Energy vs Distance';

    % Parameters
    width = 5;     % Width in inches 
    height = 3;    % Height in inches 
    alw = 1.25;    % AxesLineWidth 
    fsz = 14;      % Fontsize 
    lw = .75;        % LineWidth 
    msz = 3;       % MarkerSize 

    figure(1); 
    pos = get(gcf, 'Position');
    set(gca, 'FontSize', fsz, 'LineWidth', alw);
    x = unique(b(:,1)+1);
    len_x = size(x,1);
    y = unique(b(:,2));
    len_y = size(y,1);
    [X,Y]= meshgrid(x,y);

    freeE = real(-(0.008*310*log(b(:,3)))); 

    freeE_noInf = freeE; freeE_noInf(~isfinite(freeE_noInf))=0;
    space = round((max(freeE_noInf)-min(freeE_noInf)/4.18));
    reshaped_freeE = reshape(freeE,len_x,len_y);

    contourf(X,Y,reshaped_freeE',space);

    %Color bar
    co=colorbar;ylabel(co,'Free Energy (kJ/mol)');
    caxis([-20 10]);

    xlabel('Distance1 (nm)','FontSize',14,'FontWeight','bold');
    ylabel('Distance2 (nm)','FontSize',14,'FontWeight','bold');
