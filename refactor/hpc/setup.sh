ssh mso@graham.computecanada.ca
cd projects/def-bolker/mso
git clone https://github.com/mac-theobio/PHAC_covid.git
git clone https://github.com/bbolker/McMasterPandemic.git
module load nixpkgs/16.09 gcc/9.3.0; module load r/4.1.0
mver=`R --version | grep "R version" | awk '{print $3}' | cut -c1-3`
mkdir -p ~/R/x86_64-pc-linux-gnu-library/$mver

Rscript PHAC_covid/setup.R
