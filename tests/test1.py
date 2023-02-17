import quakeio
import ssid
import sys
from pathlib import Path
import numpy as np
from ssid import okid, srim, ExtractModes

# hwd_configs = {
#     # Case 1
#     1: {
#         # Direction (transverse->transverse)
#         "dof": {
#             "dofi": [2],
#             "dofo": [2]
#         },
#         # config dict
#         "config": {
#             "bridge":  "hayward",
#             "inputs":  [2,7,25,18],
#             "outputs": [13,15,23,20]
#         }
#     },
#     # Case 2
#     2: {
#         # Direction (longitudinal->longitudinal)
#         "dof": {
#             "dofi": [1],
#             "dofo": [1]
#         },
#         # config dict
#         "config": {
#             "bridge":  "hayward",
#             "inputs":  [3,6,17],
#             "outputs": [12,14,19]
#         }
#     },
#     # Case 3
#     3: {
#         # Direction (longitudinal->vertical)
#         "dof": {
#             "dofi": [1],
#             "dofo": [3]
#         },
#         # config dict
#         "config": {
#             "bridge":  "hayward",
#             "inputs":  [3,6,24,17],
#             "outputs": [16]
#         },
#     },
#     # Case 4
#     4: {
#         # Direction (longitudinal->longitudinal+vertical)
#         "dof": {
#             "dofi": [1],
#             "dofo": [1,3]
#         },
#         # config dict
#         "config": {
#             "bridge":  "hayward",
#             "inputs":  [3,6,24,17],
#             "outputs": [12,14,22,19,16]
#         }
#     },
#     # Case 5
#     5: {
#         # Direction (vertical->vertical)
#         "dof": {
#             "dofi": [3],
#             "dofo": [3]
#         },
#         # config dict
#         "config": {
#             "bridge":  "hayward",
#             "inputs":  [1],
#             "outputs": [16]
#         }
#     }
# }


# config = hwd_configs[1]["config"]


config = dict(
    # Transverse
    bridge  = "painter",
    # inputs  = [3],
    # outputs = [7]    
    inputs  = [17, 3, 20],
    outputs = [ 9, 7 , 4]
)


mro = 10
# mro = 20
# mro = 100
orm = 4
# orm = 10
# orm = 50
kmax = 500


# print(Path("/mnt/d/CSMIP/painter").glob("*.zip"))
print(list((Path("./")/config["bridge"]).glob("R*.zip")))



for file in (Path("./")/config["bridge"]).glob("R*.zip"):
    print(f"{file.name}")
    event = quakeio.read(file)
    try:
        inputs = [
            event.match("l", station_channel=str(chan)).accel
            for chan in config["inputs"]
        ]
        outpts = [
            event.match("l", station_channel=str(chan)).accel
            for chan in config["outputs"]
        ]
        dt = inputs[0]["time_step"]
        # A,B,C,D = ssid.srim(inputs, outpts, dt=dt, mro=100,orm=10,verbose=True)
        if True:
            # local_vars = ssid.okid.okid(inputs, outpts, dt=dt, mro=mro,orm=orm,kmax=kmax, verbose=True, debug=True)
            local_vars = ssid.srim(inputs, outpts, dt=dt, mro=mro,orm=orm, verbose=True, debug=True)
            local_vars["debug"] = False
            emacof = ssid.validation.OutputEMAC(kmax=mro,**local_vars)
            break
        A,B,C,D = ssid.okid.okid(inputs, outpts, dt=dt, mro=mro,orm=orm,kmax=kmax, verbose=True)

    except:
        raise 
        print(f"failed {file.name}")
        
    else:
        # freqdmp, modeshape, *_ = ExtractModes.ComposeModes(dt, A, B, C, D)
        # Add validation
        print(ssid.IdentifiedSystem(dt, A, B, C, D))
        # print(file, np.real(1/freqdmpSRIM[:,0]))