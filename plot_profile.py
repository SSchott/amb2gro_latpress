#!/usr/bin/env python2

import os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt

plt.rcParams['font.family'] = 'Arial'

components = ["internal","bonded","nonbonded","profile","ewald"]

def plot_profile(namd_out):
    namd_file = open(namd_out, "r")
    slab = 1
    
#    terms = {"internal":[],"bonded":[],"nonbonded":[],"profile":[]}
    terms = {component:[] for component in components}
    
    for line in namd_file.readlines():
        if line.startswith("PPROFILEBONDED:"):
            terms["bonded"].append(line.split()[2:])
        elif line.startswith("PPROFILENONBONDED:"):
            terms["nonbonded"].append(line.split()[2:])
        elif line.startswith("PPROFILEINTERNAL:"):
            terms["internal"].append(line.split()[2:])
        elif line.startswith("PRESSUREPROFILE:"):
            terms["profile"].append(line.split()[2:])
        elif "SLAB THICKNESS:" in line:
            slab = float(line.split()[-1])

    ewald_out = namd_out.replace(".namd.out","_ewald.namd.out")
    ewald = False

    if os.path.exists(ewald_out):
        ewald = True
        ewald_file = open(ewald_out, "r")
        for line in ewald_file.readlines():
            if line.startswith("PRESSUREPROFILE:"):
                terms["ewald"].append(line.split()[2:])    
    
    
    for term in terms.keys():
        print "Stats for: "+term+" in file "+namd_out
        terms[term]         = np.array(terms[term]).astype(float).reshape((len(terms[term]),len(terms[term][0])/3,3))
        terms[term+"_mean"] = terms[term].mean(axis=0)
        terms[term+"_std"]  = terms[term].std(axis=0)
    
        lat      = (terms[term][:,:,0]+terms[term][:,:,1])/2-terms[term][:,:,2]
        #lat      = (terms[term][:,:,0]+terms[term][:,:,1])/2-1
        lat_mean = lat.mean(axis=0)
        lat_std  = lat.std(axis=0)

        terms[term+"_lat"]      = lat
        terms[term+"_lat_mean"] = lat_mean
        terms[term+"_lat_std"]  = lat_std

#        z_axis = (np.arange(lat_mean.shape[0])+0.5)*slab
#        plt.fill_between(z_axis-z_axis[-1]/2,lat_mean-lat_std,lat_mean+lat_std)
#        plt.plot(z_axis-z_axis[-1]/2,lat_mean, color="r")
#        plt.savefig(term+".png")
#        plt.close()
    if ewald:
        term = "profile"
        #terms[term]         = sum([terms[t][1:] for t in components[:-2]])+terms["ewald"]
        terms[term]         = terms["profile"][1:]+terms["ewald"]
        terms[term+"_mean"] = terms[term].mean(axis=0)
        terms[term+"_std"]  = terms[term].std(axis=0)
    
        lat      = (terms[term][:,:,0]+terms[term][:,:,1])/2-terms[term][:,:,2]
        print np.mean(terms[term][:,:,2])
        print np.std(terms[term][:,:,2])
        #lat      = (terms[term][:,:,0]+terms[term][:,:,1])/2-1
        lat_mean = lat.mean(axis=0)
        lat_std  = lat.std(axis=0)

        terms[term+"_lat"]      = lat
        terms[term+"_lat_mean"] = lat_mean
        terms[term+"_lat_std"]  = lat_std

    
    return terms, slab


results = []

if len(sys.argv) ==2:
    plot_profile(sys.argv[1])
else:
    for namd_out in sys.argv[1:]:
        results.append(plot_profile(namd_out))

    means_mean         = {component:[] for component in components}
    means_stderr       = {component:[] for component in components}
    means_stderr_prop  = {component:[] for component in components}
    means              = {component:[] for component in components}
    
    for component in components:
        
        means_mean[component]         = np.array([result[0][term] for result in results for term in result[0] if term == component+"_lat_mean"]).mean(axis=0)
        means_stderr[component]       = np.array([result[0][term] for result in results for term in result[0] if term == component+"_lat_mean"]).std(axis=0)/np.sqrt(len(results))
        means_stderr_prop[component]  = np.sqrt(np.sum(np.array([result[0][term]**2 for result in results for term in result[0] if term == component+"_lat_std"]), axis=0))/(len(results)*np.sqrt(len(results)))
        means[component]              = np.array([result[0][term] for result in results for term in result[0] if term == component+"_lat_mean"])

        z_axis = (np.arange(results[0][0][component+"_lat_mean"].shape[0])+0.5)*results[0][1]
        plot_axis = z_axis-z_axis[-1]/2
        
    
        plt.fill_between(plot_axis, means_mean[component]-means_stderr[component],means_mean[component]+means_stderr[component])
        plt.plot(plot_axis,means_mean[component], color="r")
        plt.xlim((min(plot_axis),max(plot_axis)))
        plt.xlabel("z-coordinate [$\AA$]")
        plt.ylabel("Pressure [bar]")
        plt.savefig(component+"_mean.png")
        plt.close()
    
        plt.fill_between(plot_axis, means_mean[component]-means_stderr_prop[component],means_mean[component]+means_stderr_prop[component])
        plt.plot(plot_axis,means_mean[component], color="r")
        plt.xlim((min(plot_axis),max(plot_axis)))
        plt.xlabel("z-coordinate [$\AA$]")
        plt.ylabel("Pressure [bar]")
        plt.savefig(component+"_mean_prop.png")
        plt.close()
    
    
        np.savetxt(component+".dat"      , np.c_[plot_axis,means_mean[component],means_stderr[component],means_stderr_prop[component]])
        np.savetxt(component+"_means.dat", np.c_[plot_axis,means[component].T])
