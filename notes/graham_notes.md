A quick update on attempting to run calibrations on Graham:
As expected, our standard "pure" Nelder-Mead operation is running normally. This completes a test calibration in 17 minutes, and all the metrics look normal. A calibration using deOptim with only one core allocated doesn't complete in 20 minutes (I can try rerunning it for longer), but has otherwise normal-looking metrics. A calibration using deOptim with 4 cores (doesn't even complete in 40 minutes, and has an extremely low CPU usage indicating that it is not being used fully. Note that I accidentally over-allocated memory for the Nelder-Mead and deOptim 1-core runs, resulting in low memory efficiency.

The calibrate call is as follows (for the 4-core deOptim):
ont_cal1 <- calibrate(
  data = ont_all_sub,
  base_params = params,
  opt_pars = opt_pars,
  time_args = list(params_timevar = timevar_df),
  sim_args = list(step_args = list()),
  use_DEoptim = TRUE,
  DE_cores = 4
)

(Pure Nelder-Mead)
[mso@gra-login3 scratch]$ seff 49577068
perl: warning: Setting locale failed.
perl: warning: Please check that your locale settings:
        LANGUAGE = (unset),
        LC_ALL = (unset),
        LANG = "C.UTF-8"
    are supported and installed on your system.
perl: warning: Falling back to the standard locale ("C").
Job ID: 49577068
Cluster: graham
User/Group: mso/mso
State: COMPLETED (exit code 0)
Cores: 1
CPU Utilized: 00:16:53
CPU Efficiency: 99.31% of 00:17:00 core-walltime
Job Wall-clock time: 00:17:00
Memory Utilized: 326.91 MB
Memory Efficiency: 7.98% of 4.00 GB

(deOptim 1-core)
[mso@gra-login3 scratch]$ seff 49577069
Job ID: 49577069
Cluster: graham
User/Group: mso/mso
State: TIMEOUT (exit code 0)
Cores: 1
CPU Utilized: 00:19:52
CPU Efficiency: 98.27% of 00:20:13 core-walltime
Job Wall-clock time: 00:20:13
Memory Utilized: 354.03 MB
Memory Efficiency: 8.64% of 4.00 GB

(deOptim 4-core)
[mso@gra-login1 scratch]$ seff 49578001
perl: warning: Setting locale failed.
perl: warning: Please check that your locale settings:
        LANGUAGE = (unset),
        LC_ALL = (unset),
        LANG = "C.UTF-8"
    are supported and installed on your system.
perl: warning: Falling back to the standard locale ("C").
Job ID: 49578001
Cluster: graham
User/Group: mso/mso
State: TIMEOUT (exit code 0)
Nodes: 1
Cores per node: 4
CPU Utilized: 00:07:32
**CPU Efficiency: 6.19% of 02:01:40 core-walltime**
Job Wall-clock time: 00:30:25
Memory Utilized: 327.13 MB
Memory Efficiency: 15.97% of 2.00 GB
