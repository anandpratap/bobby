import sys, os
import numpy as np
import matplotlib.pyplot as plt

sys.path.append(os.path.join(os.path.dirname(__file__), "../"))
from src.bobby import Bobby

b = Bobby()
b.run()
