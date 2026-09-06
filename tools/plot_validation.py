"""Export scientific plots and CSVs from the two actual complete case runs."""
from pathlib import Path
import argparse
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'tests'))
from compare_full_case import acoustic,CASE

def main():
    p=argparse.ArgumentParser();p.add_argument('baseline',type=Path);p.add_argument('candidate',type=Path)
    args=p.parse_args();output=ROOT/'examples/results';output.mkdir(parents=True,exist_ok=True)
    a,_,_=acoustic(args.baseline/f'{CASE}_1.out');b,_,_=acoustic(args.candidate/f'{CASE}_1.out')
    sa,header,_=acoustic(args.baseline/f'{CASE}_2.out');sb,_,_=acoustic(args.candidate/f'{CASE}_2.out')
    frequencies=np.array([float(x.split('Freq')[-1]) for x in header[1:35]])
    mean_a=10*np.log10(np.mean(10**(sa[:,1:]/10),axis=0)).reshape(2,34)
    mean_b=10*np.log10(np.mean(10**(sb[:,1:]/10),axis=0)).reshape(2,34)
    np.savetxt(output/'official_oaspl.csv',np.column_stack((a,b[:,1:])),delimiter=',',comments='',
        header='time_s,fortran_obs1_dB,fortran_obs2_dB,cpp_obs1_dB,cpp_obs2_dB',fmt='%.9g')
    np.savetxt(output/'official_mean_spectrum.csv',np.column_stack((frequencies,mean_a.T,mean_b.T)),delimiter=',',comments='',
        header='frequency_Hz,fortran_obs1_dB,fortran_obs2_dB,cpp_obs1_dB,cpp_obs2_dB',fmt='%.12g')
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig,axes=plt.subplots(1,2,figsize=(11,4.2),layout='constrained')
    colors=['#006b92','#cf641f']
    for j,color in enumerate(colors):
        axes[0].plot(a[:,0],a[:,j+1],color=color,lw=2,label=f'Observer {j+1}, Fortran')
        axes[0].plot(b[::10,0],b[::10,j+1],color=color,ls='none',marker='o',ms=4,mfc='white',label=f'Observer {j+1}, C++')
        axes[1].semilogx(frequencies,mean_a[j],color=color,lw=2)
        axes[1].semilogx(frequencies,mean_b[j],color=color,ls='none',marker='o',ms=4,mfc='white')
    axes[0].set(xlabel='Time (s)',ylabel='OASPL (dB re 20 uPa)',title='Full 20 s simulation',xlim=(0,20))
    axes[0].legend(fontsize=8,ncols=2)
    axes[1].set(xlabel='Frequency (Hz)',ylabel='Band SPL (dB re 20 uPa)',title='Energy-averaged spectrum, 0–20 s',xlim=(10,20000))
    for ax in axes:ax.grid(True,alpha=.2)
    fig.suptitle('IEA_LB_RWT-AeroAcoustics | GCC 16.2 + Intel oneMKL',fontsize=13)
    fig.savefig(ROOT/'docs/full-case-validation.png',dpi=160)
    print('Exported actual case comparison plot and CSVs')

if __name__=='__main__':main()
