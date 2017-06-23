# L1HATS
Git repository for L1 HATS tutorial

### Setup CMSSW project area
<pre><code>
setenv SCRAM_ARCH slc6_amd64_gcc530
cmsrel CMSSW_9_2_2
cd CMSSW_9_2_2/src
cmsenv
</code></pre>

### L1Menu DPG package for menu-making
<pre><code>
git clone https://github.com/cms-l1-dpg/L1Menu.git L1TriggerDPG/L1Menu
cd L1TriggerDPG/L1Menu
scramv1 b -j 9
cd macros
make -j 9
</code></pre>

### HATS L1 Tutorial
<pre><code>
cd $CMSSSW_BASE/src
git clone git@github.com:apana/L1HATS.git
cd L1HATS/L1HATSexercise/macros
make
</code></pre>

### Run the example and plot the rates
<pre><code>
rateExample -n 1000000 -b 1549 -o Rates.root -l run297100_ZeroBias_noRECO.list
root Rates.root
</code></pre>
