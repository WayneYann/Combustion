
from hopper import HopperPBS
from serial import Serial

by_name = {
    'hopper': HopperPBS,
    'edison': HopperPBS,
    'serial': Serial,
}
