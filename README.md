# WIMPyCCD

These scripts calculate the nuclear recoil interaction rate between dark matter particles (of a certain mass and DM-n cross-section) and a detector target of a given material of mass number A. This is done analytically, by considering the movement of Earth and Sun in the galaxy.

The ionization inefficiency effects of different models (DAMICM 2025/SuperCDMS measurement, Lindhard model, etc.) are implemented to transform the nuclear recoil energy to the ionized electron's energy (i.e., visible energy). This has significant effects on low-energy events. The output is a differential energy spectrum of ionized electrons.

The energy spectrum of ionized electrons is calculated for each dark matter mass and dm-n cross-section at a given detector material, threshold energy to extract the projected upper limit of permitted DM-n cross-section & DM mass phase space at 90% confidence level if no event is observed during certain mass-time exposure (assume no background for now).

To dos: Pixelization, detector response, reconstruction efficiency, and background effects
(which is what I hoped to talk about. Since now I know how to implement many of them via wader, and all I need is to generate event-based simulations of detector hits (with or without Geant4?)
Lewin and Smith: https://www.sciencedirect.com/science/article/pii/S0927650596000473

Ben Loer's Thesis: https://borex.lngs.infn.it/Thesis/B.M.Loer_PhD_Thesis.pdf
