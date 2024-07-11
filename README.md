

### Author Contact Information

russell.jasper@unibe.ch
rjjasper@ucalgary.ca

### Usage and license information

If you use or are inspired by code in this repository please cite the following work or contact me about how to cite. 

Jasper & Yeaman 2020 Local adaptation can cause both peaks and troughs in nucleotide diversity within populations. bioRxiv doi: https://doi.org/10.1101/2020.06.03.132662

---

# Local adaptation can cause both peaks and troughs in nucleotide diversity within populations

I used Nemo and R to perform forward in time simulations in order to investigate how local adaptation shapes patterns in nucleotide diversity and Fst at neutral sites linked to a causal locus. I explored how migration, mutation, degrees of polygenicity & redundancy, and population structure affected variation at linked neutral loci.

---

# Environment

Nemo version 2.3.46

R version 1.4.1106

# Usage

### 01_init_maker

This code outputs both .init and .slurm files to run Nemo simulations for multiple different parameter sets on the slurm server of your choice.

Edit the "init_template.txt" file with the simulation parameters relevant to you. This will be the static foundation for all simulation reps, we will add the variable parameters next. 

In R:

  -Read in the init_template.txt
  
  -Input the values for the variable parameters you wish to iterate over
  
  -Input relevant slurm parameters, eg path to Nemo program and time, mem, etc
  

Run loop --> collect .init files and .slurm files.

Place files in the same directory on slurm server. Submit the .slurm file to run the corresponding .init file for each particular parameter set.

### 02A_calculate_pi / 02B_calculate_fst

This code takes the resulting .txt (or .freq) outputs from Nemo and calculates the mean pi (or Fst) for each locus in a patch (or metapopulation) for all parameter sets.

In R, point to the directory containing the relevant raw Nemo results files and input a string for the DATE for naming the output files.

Run code --> collect .csv files containing pi or Fst.

### 03_plot_pi_fst

This code takes the previous results and plots with ggplot.

Figure 1 (see: Jasper and Yeaman 2020 Fig. 1A)
On the y-axis I have used the slope of a metric against the distance from the adaptive locus, ie, a positive slope suggests a depression in the metric nearby the adaptive locus and a negative slope suggests an inflation in the metric nearby the adaptive locus. The x-axis and the colour are the migration rate in log space, and the variance in the fitness function, respectively. Adjust as you see fit.

Figure 2 (see: Jasper and Yeaman 2020 Fig. 1B)
On the y-axis I have plotted the metric at two loci, one 1/1000 of a cM away from the adaptive locus and one 10 cM away from the adaptive locus. Other aspects are as described above.

Figure 3
On the y-axis it is simply the per locus mean of the metric across the entire chromosome. Other aspects are as described above.


# Repo Table of Contents

## Data

### example_init_template
An example init_template.txt to generate .init and .slurm files for running simulations.

## Code
As described above.

### 01_init_maker

### 02A_calculate_pi / 02B_calculate_fst

### 03_plot_pi_fst



