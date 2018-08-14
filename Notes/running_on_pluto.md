# Tips for running things on Pluto (the workstation)

## Variables
  - `$PID`: Process ID. Should be an integer of more than 4 or so digits. Can be found using the `jobs` command.
  - `$JID`: Job ID. Should be an integer of single digits. Can be found using the `jobs` command.
  - `$UNAME`: Your username (should be your JHID).
  - `$PROC_NAME`: Name of a process.

## General
  - use the `nohup` command.
    - usage: `nohup command`
    - `nohup` stands for "no hangup", meaning if the terminal/ssh session is killed, the process does not get killed
    - if you're using any shell other than `fish`, and you forget to use `nohup`, you can retroactively mimic it by using `disown $JID`
  - To pause a job, use `kill -TSTP $PID` and continue it with `kill -CONT $PID`
  - to see what a particular terminal is running, type `jobs`
    - jobs in a terminal can be killed with `jobs -9 %$JID`.
  - to see what you're running, if you cannot use the `jobs` command, type `ps -u $UNAME`. This will output every process your username is currently in ownership of
    - can whittle this down with `pgrep $PROC_NAME -u $UNAME | xargs --no-run-if-empty ps fp`.
    - for example, if I wanted to see what Gromacs processes I am running, I would type `pgrep gmx -u mvarga3 | xargs --no-run-if-empty ps fp`
  - If you want to log everything Gromacs outputs, type `> $LOG 2>&1` at the end of the MD run command,
    ```bash
    # This will log everything mdrun outputs
    # 1 stands for all standard output, 2 stands for standard error output
    # so 2>&1 outputs all standard output to a file, while concatenating all standard error to the same file
    gmx mdrun -v -deffnm md > md.log 2>&1
    ```
  - When attempting to run VMD over an ssh session while using macOS, if you get an error along the lines of a bad integer value for XRequest not allowing it to create an OpenGL window, follow these steps:
    1. Quit the ssh session
    2. Type `defaults write org.macosforge.xquartz.X11 enable_iglx -bool true` 
    3. Type `reset` or kill the terminal and open a new one. You should now be able to open a VMD window from the remote machine. Viewing trajectories in it does not work well, though.

### Process Priority
  - Niceness dictates what process a CPU will prioritize, on a scale [-20,20), with lower niceness being prioritized over higher niceness. This means that a process with niceness 0 will be prioritized, in terms of CPU time, over a process with niceness 1 (niceness can be seen with commands such as `ps` and `top`).
  - To set niceness at the beginning of a process, use `nice -n $NICE $CMD`, where `$NICE` is the niceness and `$CMD` is whatever command you want to run.
  - To set the niceness of a currently running process, use `sudo renice $NICE -p $PID`.

## Gromacs with Plumed
  - Plumed needs to be compiled with MPI enabled (`--enable-mpi`)
  - Plumed needs to have at least one environment variable set:
    - `set PLUMED_KERNEL=$PATH`, for whatever the Plumed kernel path is (default should be `/usr/local/lib/libmpi-plumed_2.3Kernel.so`)
  - use `gmx-plumed2.3` instead of `gmx`
  - When running `gmx-plumed2.3 mdrun`, include the flags `-ntmpi 1 -ntomp 16`, otherwise it will not run
