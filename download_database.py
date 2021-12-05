#starting importing some useful stuff
import os

from pypath.share import settings
settings.setup(progressbars = True)

#setting the cache directory
os.makedirs('tmpcache_network')
settings.setup(cachedir = 'tmpcache_network') # ==> to use when network.Network will be fuly implemented


#importing legacy, which is the 'old' version of pypath, the only one (for now), with the graph object implemented
from pypath.legacy import main as legacy
#importing the new class used by pypath, useful for the dataframe visualization
#but with noreferences (maybe a bug to report?). Here graph is not implemented
from pypath.core import network
from pypath.resources import network as netres

#initialization of the network object
pa = network.Network()

#this is used to build the network with the new class Network
#unfortunately the resulting network cannot be visualized with igraph and there re no reference
for database in netres.data_formats.omnipath.keys():
    if (database=="netpath"):
        continue
    try:
        lst={database: netres.data_formats.omnipath[database]}
        print(database, " : ", netres.data_formats.omnipath[database])
        pa.init_network(netres.data_formats.omnipath[database])
    except Exception as inst:
        print("Error for "+database)

#since build the network is super expansive in terms of time, once built the network it is possible to save it
#as pickle file and load it in less time

## save the network:
pa.save_to_pickle(pickle_file='mynetwork.pickle')

print('END OF THE DOWNLOAD/UPDATE OF THE NETWORK/n')



os.makedirs('tmpcache_legacy')
settings.setup(cachedir = 'tmpcache_legacy') # ==> actual cache folder to use with legacy.main

#initialization of the 'old' PyPath object
pw_legacy = legacy.PyPath()

#USING LEGACY

#DO NOT EVALUATE UNLESS YOU WANT TO UPDATE THE DATABASES

#instead of loading the databases using steps, I prefer loading the activity flow networks with literature references
#(you can look at the possible datasets at: https://workflows.omnipathdb.org/pypath_guide.html#network-resources)

for database in legacy.data_formats.omnipath.keys():
    try:
        print(database, " : ", legacy.data_formats.omnipath[database])
        lst={database: legacy.data_formats.omnipath[database]}
        pw_legacy.init_network(lst)
    except Exception as inst:
        print("Error for "+database)

for database in legacy.data_formats.ligand_receptor.keys():
    try:
        print(database, " : ", legacy.data_formats.ligand_receptor[database])
        lst={database: legacy.data_formats.ligand_receptor[database]}
        pw_legacy.init_network(lst)
    except Exception as inst:
        print("Error for "+database)

for database in legacy.data_formats.tf_mirna.keys():
    try:
        print(database, " : ", legacy.data_formats.tf_mirna[database])
        lst={database: legacy.data_formats.tf_mirna[database]}
        pw_legacy.init_network(lst)
    except Exception as inst:
        print("Error for "+database)


#it still requires lots of time so better save it

#save the network

pw_legacy.save_to_pickle('mylegacy.pickle')

import sys
sys.exit("END OF THE SCRIPT")