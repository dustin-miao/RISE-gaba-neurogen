In this directory the models used for Fig. 4 of the manuscript 
"NKCC-1 mediated Cl- uptake in immature CA3 pyramidal neurons is 
sufficient to compensate phasic GABAergic input" by Kolbaev et al.
submitted to Scientific Reports

To use this model please unpack in folder, and use mknrndll to 
compile the used .mod files.

For the simulation used in Fig. 4B:
  Start Neuron by clicking Cell_1_SciRep_ShrinkCorr.hoc
  Load session "Determine_tau_NKCC1_rig.ses"
  In the cldif_CA3_NKCC1_HCO3 window you can adjust the Values for
  the NKCC1 mediated Cl- transport.
  Most relevant are the parameters: 
    tau_NKCC1 - adjust this value to fir the measured Cl- trajectory
    cli_Start - adjust to the desired Starting Value, 9.1 mM is the exp. determined value

For the simulation used in Fig. 4C:
  Start Neuron by clicking Cell_1_SciRep_ShrinkCorr.hoc
  Load session "Determine_cl-Flux_w-o_NKCC1_rig.ses"
  In the cldif_CA3_NKCC1_HCO3 window you can adjust the Values for
  the NKCC1 mediated Cl- transport.
  Most relevant are the parameters: 
    tau_NKCC1 - adjusted to 9.9*10(19) i.e. pratically absent 
    tau_passive - adjusted to 9.9*10(19) i.e. pratically absent 
                  (note that in this paper we evaluate the Cl- fluxes 
                  via the tonic GABA A receptors implemented in Cell1)
    cli_Start - adjust to the desired starting Value, here 17 mM as this is the 
                steady state of the recordings shown in Fig. 4C

For the simulation used in Fig. 4E:
  Start Neuron by clicking start_Phasic_GABA_activity_only_soma_Backregul.hoc

For the simulation used in Fig. 4F:
  Start Neuron by clicking start_Phasic_GABA_activity_only_soma_Div_gGABA.hoc
  The different conductances  of the phasic currents are sequentially
  played in the main file "Phasic_GABA_activity_only_soma_Div_gGABA.hoc" and
  stored in .asc files

For the simulation used in Fig. 4G:
  Start Neuron by clicking start_Phasic_GABA_activity_only_soma_Div_Freq.hoc
  The different conductances  of the phasic currents are sequentially
  played in the main file "Phasic_GABA_activity_only_soma_Div_gGABA.hoc" and
  stored in .asc files Start-Frequency and number of frequency staps are defined
  in line 15 and 16 of the hoc file. For the manuscript thes are set manually
  by defining Freq-Staps = 1 and Start_Freq to 1,2,5,10,20, and 50
  
For the simulation used in Fig. 4H:
  Start Neuron by clicking "start_Block_Tonic_Cl-currents.hoc"
 
For the simulation used in Fig. 4I:
  Start Neuron by clicking "start_Add_Tonic_Cl-currents.hoc"
  The amplitude of the added tonic current mut be set in the Valuable Add_Tonic
  in the second line of the code in "Add_tonic_Cl-current.hoc
 


