#  Weekly Calibrations
This notebook contains the necessary code to set up the weekly regression testing framework, allowing us to check if the internals have been significantly broken by a given change. Please change my credentials to whatever yours are.

### Working on Graham 
I guess any similar Compute Canada cluster should do, but we're going to use SHARCNET/Graham.
```
ssh mso@graham.computecanada.ca
```

Create your ssh key on Graham, you will use this to ssh into GitHub servers.
```
ssh-keygen -t rsa -b 4096 -C your_email@example.com
```

Copy-paste your RSA public key. After you have it, copy-paste it into GitHub following [these](https://docs.github.com/en/github/authenticating-to-github/connecting-to-github-with-ssh/adding-a-new-ssh-key-to-your-github-account) steps.
```
cd ~/.ssh
vim id_rsa.pub # or similar, wherever you put your key!
```

Test your ssh connection to GitHub.
```
ssh -T git@github.com
```

You will now clone and/or install McMasterPandemic. Follow the instructions in the [README](https://github.com/mac-theobio/McMasterPandemic) to get set up with MacPan if you haven't done so already. 
```
cd /home/mso/projects/def-bolker/mso
git clone ssh://git@github.com/mac-theobio/McMasterPandemic.git
```

Let's leave Graham and move to another server.
```
exit
```
### Working on a scheduler server

Again, we ssh into a server with the cron scheduler. Graham doesn't offer cron, so we have to ssh Graham from another server (in this case, the mathematics and statistics server.)
```
ssh som5@ms.mcmaster.ca # change to math server
ssh-keygen -t rsa -b 4096 -C your_email@example.com # create authentication key
ssh-copy-id mso@graham.computecanada.ca # instead of copy-pasting your key to GitHub, we can just use a command to move it to Graham.
ssh mso@graham.computecanada.ca 'cd /home/mso/projects/def-bolker/mso/McMasterPandemic/misc && touch test_file.txt' # execute a test ssh command
```

Then, we set up our job:
```
service crond status # check cron status 
crontab -e
```

This opens a vim session. Paste this in here, modify the cd command to match your directory structure. This cron job executes at 12:00 AM every Monday.
```
#!/bin/bash
PATH=/usr/bin
MAILTO=your_email_here@example.com

0 0 * * 1 ssh mso@graham.computecanada.ca "cd /home/mso/projects/def-bolker/mso/McMasterPandemic/misc && git add --all && git commit -am "calibration sync" && git pull && sbatch calibration_sbatch.sh"
```

And you should be set up now!

### How it works
The actual scripts might change a bit from these commands, but the overall explanation should remain the same.

Let's break-down the crontab ssh command.
- ```ssh mso@graham.computecanada.ca```: Run this command on Graham.
- `cd /home/mso/projects/def-bolker/mso/McMasterPandemic/misc`: Change working directory to the MacPan repo.
- `sbatch calibration_sbatch.sh`: Graham syntax for scheduling a script to work in the future.

Then, let's look at the `calibration_sbatch.sh`.

- These are typical comments needed to specify the job.
  - ``` #!/bin/bash
      #SBATCH --time=00:25:00
      #SBATCH --account=def-bolker
      #SBATCH --job-name=ont_cal_NelderMead
      #SBATCH --cpus-per-task=1
      #SBATCH --mem-per-cpu=400M
    ```
-  `module load r/4.1.0`: Load R.
-  `loc=$(pwd); ssh gra-login1 "ssh gra-login1 "cd $loc && git add --all && git commit -am 'calibration sync' && git pull"`: ssh the login node from the compute node (on which the job is run). The compute node has no internet access, but can ssh the login node (which does). This way, we can sync the local repo with origin.
-  `Rscript ont_cal_NelderMead.R`: Our calibration script.
-  `loc=$(pwd); ssh gra-login1 "cd $loc && bash git_push.sh`: Do a similar trick, but to push our outputs.
