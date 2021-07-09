ssh mso@graham.computecanada.ca # log into Graham
cd /home/mso/projects/def-bolker/mso
ssh-keygen -t rsa -b 4096 -C somatthewc@gmail.com
cd ~
# find ssh key and add to GitHub https://docs.github.com/en/github/authenticating-to-github/connecting-to-github-with-ssh/generating-a-new-ssh-key-and-adding-it-to-the-ssh-agent
git clone ssh://git@github.com/mac-theobio/McMasterPandemic.git
# I had already installed the stuff necessary for McMasterPandemic, but presumably you would do that here.
# create sbatch work script (see misc/calibration_sbatch.sh), this runs on compute clusters
# create submission script that waits for sbatch script to finish, see misc/weekly_calibrations.sh



ssh som5@ms.mcmaster.ca # change servers to math server
ssh-keygen -t rsa -b 4096 -C somatthewc@gmail.com
ssh-copy-id mso@graham.computecanada.ca
ssh mso@graham.computecanada.ca 'cd /home/mso/projects/def-bolker/mso/McMasterPandemic/misc && bash push_test.sh'
service crond status # check cron status 
crontab -e

# put these lines in it for testing
# ----------------------------------------------------------------------------------------------------------------
#!/bin/bash
PATH=/usr/bin
MAILTO=somatthewc@gmail.com

* * * * * ssh mso@graham.computecanada.ca 'cd /home/mso/projects/def-bolker/mso/McMasterPandemic/misc && bash push_test.sh'
# ----------------------------------------------------------------------------------------------------------------

# this should be the real crontab
# ----------------------------------------------------------------------------------------------------------------
#!/bin/bash
PATH=/usr/bin
MAILTO=somatthewc@gmail.com

* * * * * ssh mso@graham.computecanada.ca 'cd /home/mso/projects/def-bolker/mso/McMasterPandemic/misc && bash weekly_calibration.sh'
# ----------------------------------------------------------------------------------------------------------------
