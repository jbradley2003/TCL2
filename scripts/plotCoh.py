import os
import pickle
import matplotlib.pyplot as plt
import numpy
from functions import *

T = 'Times (millisecond)'
R = 'Coherence'

# with open(f'{os.getcwd()}/tempo_bL30_nB1000_c10_TCL2_filtered_geminal_no_methyl_no_bath.pkl', 'rb') as f:
#         tempoNB = pickle.load(f)

# with open(f'{os.getcwd()}/tempo_bL30_nB1000_c10_TCL2_filtered_geminal_no_bath.pkl', 'rb') as f:
#         tempoJNB = pickle.load(f)

with open(f'{os.getcwd()}/tempo_bL30_nB100_c10_TCL2_filtered_geminal_no_methyl.pkl', 'rb') as f:
        tempo = pickle.load(f)

with open(f'{os.getcwd()}/tempo_bL30_nB100_c10_TCL2_filtered_geminal.pkl', 'rb') as f:
        tempoJ = pickle.load(f)


# plt.plot(tempo[T], tempo[R])
# plt.plot(tempo[T], tempoJ[R])
plt.xlabel('Time (ms)')
plt.ylabel(r'$\rho_e^{01}$')
plt.title('Bath Dynamics / bL = 30, nB = 25')
# plt.title('Single Molecule Dynamics (noBath) / bL = 30, nB = 1000')
plt.tick_params(direction='out', length=6, width=1)
a, b = fit_exp(tempo[R],tempo[T]*1e3)[1][0], fit_exp(tempoJ[R],tempoJ[T]*1e3)[1][0]
plt.plot(tempo[T]*1e3, tempo[R], label=r'$\gamma_{kl} = b_{kl}$, $T_2=$'f'{np.round(a,3)}'r'$\mu s$')
plt.plot(tempoJ[T]*1e3, tempoJ[R],label=r'$\gamma_{kl} = b_{kl} - 2J_{kl}$, $T_2=$'f'{np.round(b,3)}'r'$\mu s$')
# plt.show()

# plt.plot(tempo[T]*1e3, tempoNB[R], label=r'$\gamma_{kl} = b_{kl}$')
# plt.plot(tempo[T]*1e3, tempoJNB[R],label=r'$\gamma_{kl} = b_{kl} - 2J_{kl}$')

plt.legend()
plt.legend(frameon=False, loc='upper right', handlelength=3)
plt.show()

