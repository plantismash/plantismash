#!/usr/bin/env python
import sys
import os
from math import log
from random import random
from numpy import array
from hmms import *
import pickle, copy
from antismash import utils
from operator import itemgetter

# link to where the frequencies are
PATH_TO_FREQUENCIES = utils.get_full_path(__file__, "")
### Read in Files !
out4 = open(PATH_TO_FREQUENCIES + 'TP_arr_A.pkl','r')
out5 = open(PATH_TO_FREQUENCIES + 'NewTS_all_B_reduced_6filter.pkl','r')
out6 = open(PATH_TO_FREQUENCIES + 'SP_arr.pkl','r')
out7 = open(PATH_TO_FREQUENCIES + 'NewTS_all_B_index.pkl','r')

A=pickle.load(out4)
B=pickle.load(out5)
start_p=pickle.load(out6)
index_dict=pickle.load(out7)

out4.close()
out5.close()
out6.close()
out7.close()

repeats = ['PF07721', 'PF05593', 'PF07719', 'PF00515', 'PF00132', 'PF03130', 'PF01839', 'PF01816', 'PF07720', 'PF00400', 'PF05594', 'PF07661', 'PF02985', 'PF06049', 'PF08238', 'PF06696', 'PF00353', 'PF02412', 'PF00023', 'PF02071', 'PF03991',     'PF01469', 'PF07676', 'PF00514', 'PF00904', 'PF07634', 'PF02370', 'PF03335', 'PF01851', 'PF04728', 'PF06715', 'PF03373', 'PF04680', 'PF00805', 'PF04508', 'PF07918', 'PF01535', 'PF01011', 'PF05017', 'PF06671', 'PF00818', 'PF03406', 'PF00399',     'PF09373', 'PF01744', 'PF01436', 'PF01239', 'PF05906', 'PF03729', 'PF00404', 'PF04022', 'PF02363', 'PF02524', 'PF07981', 'PF02095', 'PF00414', 'PF00560', 'PF05001', 'PF02162', 'PF01473', 'PF05465', 'PF02493', 'PF03578', 'PF08043', 'PF06392',     'PF07142', 'PF08309', 'PF02184']

# --- HMM run
def do_HMM(observations,A,B,index_dict,start_p,graph=False,report=False,outs=None,method='HMM',path=os.getcwd()+"/",name='bla',FilterRepeats=0,BasesPerLine = 500000):

  '''
  This is an interface between HMM module and users data. It takes in Python array of Pfam IDs ordered by their genomic location, and returns an array
  of gene cluster probabilities. Optionally, it can also write a text report and plot a lousy probabilites versus pfam position graph.

  Inputs:
    Mandatory:
    - observations: Python array of Pfam IDs (example: ['PF00001','PF00550','PF00501',...])
    - A: array of transition probabilities (read HMM module documentation for more info)
    - B: array of emission probabilites (read HMM module documentation for more info)
    - index_dictio: dictionary of Pfam indeces for positions in A
    - start_p: vector of start probabilities
    Optional:
    - graph(0/1): plot graph (0 by default)
    - report(0/1): write report (0 by default)
    - outs: Pyton array of arrays with information about each Pfam you would like to write in output,
          a probability will be append to each of the sub-arrays, tab-delimited, and written in the output file
    - method: HMM for Forward-Backward algorithm, viterbi for Viterbi alogirthm (HMM by default)
    - path: path to input and output files
    - name of the output files (graph and/or report)
                - filter_repeats: pass to ploter if you're filtering repeats
                - BasesPerLine: if plotting, the number of basepairs per line
  '''

  index_dictio = copy.deepcopy(index_dict)

  reals = None
  if method=='viterbi':
    print >>sys.stderr, '\tPredicting Clusters...'
    PREDICTED = list((array(viterbi(A,B,start_p,observations,index_dictio))-1)*(-1))
    print >>sys.stderr, '\tCluster predicted!\n'

  elif method == 'HMM':
    #print >>sys.stderr, '\tPredicting Clusters...'
    PREDICTED = list(forward_backward(A,B,start_p,observations,index_dict,method='BF'))
    #print >>sys.stderr, '\tCluster predicted!\n'

  else:
    print >>sys.stderr, "\nMethod: '%s' unknown. Please choose HMM or viterbi\n" % method
    return None


  if graph == True:
    Data = read_input(outs,PREDICTED,filter_repeats=FilterRepeats)
    X,Y = {},{}
    for n,r in enumerate(Data):
      x,y = data_prepare(Data[r],LineLength=BasesPerLine)
      X[n] = x
      Y[n] = y

    ploter(X,Y,name=name,path=path,LineLength=BasesPerLine)

  if report == True:
    output = open(path+name+'.out','w')
    if len(outs) == len(PREDICTED):
      for n,i in enumerate(PREDICTED):
        outs[n].append(PREDICTED[n])
        lines = '\t'.join([str(i) for i in outs[n]]) + '\n'
        output.write(lines)
    else:
      print >>sys.stderr, 'OutputWriteError: PREDICTED array length does not match wigh Genome class!'
      return None
    output.close()

  #print '%s Done!' % name
  return PREDICTED


# --- plotter
def ploter(Xaxis,Yaxis,LineLength = 500000,path=os.getcwd()+"/",name='bla'):
  '''
  Graph plotting module. Very basic, so feel free to modify it.
  '''
  try:
    from matplotlib.figure import Figure
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.ticker import MultipleLocator
  except ImportError:
    print >>sys.stderr, '\nMatplotlib is not installed! Either install it, or deselect graph option.\n'
    return 1

  fontsize = 8
  Lines = sum([len(Xaxis[i]) for i in Xaxis])

  fig = Figure(figsize=(Lines,16),dpi=100)
  fig.subplots_adjust(hspace=0.4)
  Colors = ['r','g','b','c','m','y','k']*10

  C = 0
  for chain in Xaxis:
    X = Xaxis[chain]
    Y = Yaxis[chain]

    for x in xrange(len(X)):
      C += 1

      ax = fig.add_subplot(Lines,1,C)
      if C == 1:
        ax.set_title('Biosynthetic Gene Cluster Probabilities')

      ax.plot(X[x],Y[x],'%s-' % Colors[chain],linewidth = 2)

      ax.set_xlim(x*LineLength,(x+1)*LineLength)
      ax.set_ylim(0.,1.)
      ax.set_xticklabels(range(x*LineLength/1000,(x+1)*LineLength/1000)[::100] + [(x+1)*LineLength/1000])
      ax.set_yticklabels(())
      ax.xaxis.set_minor_locator(MultipleLocator(10000))

      for tick in ax.xaxis.get_major_ticks():
        tick.label1.set_fontsize(fontsize)
      for tick in ax.yaxis.get_major_ticks():
        tick.label1.set_fontsize(fontsize)
        ax.set_xlabel('Genomic position (kbp)')

        canvas = FigureCanvasAgg(fig)
        canvas.print_figure(path + os.sep + name + '.cluster_prediction.png', dpi=100)

# --- prepare input list for plotting
def data_prepare(A,LineLength = 500000):

  X,Y = [],[]
  XXX,YYY = [],[]
  L = LineLength
  for x in xrange(0,len(A)):
    if A[x][0] < L:
      X.append( (A[x][0] + A[x][1])/2 )
      Y.append( A[x][3] )
    else:
      X.append( L )
      Y.append( A[x][3] )
      L += LineLength
      XXX.append(X)
      YYY.append(Y)
      X,Y = [],[]

  if len(X) != 0:
    XXX.append(X)
    YYY.append(Y)

  return (XXX,YYY)

# --- read input list
def read_input(outs,predictions,filter_repeats=0):

  global repeats

  Info = {}
  x = 0
  for d in outs:

    # change the next 4 lines according to Predicter's output file
    chain,geneID = 'chain1',None
    gStart,gStop = None,None
    PfamID,prob = d[0],predictions[x]
    pStart,pStop = int(d[2]),int(d[3])

    if filter_repeats==1:
      if d[0] in repeats: continue

    if chain not in Info: Info[chain] = [[pStart,pStop,PfamID,prob]]
    else: Info[chain].append([pStart,pStop,PfamID,prob])

    x += 1

  Info1 = {}
  for chain in Info.keys():
    Info1[chain] = sorted(Info[chain],key=itemgetter(0))

  Info = None
  return Info1

# --- filter repeats
def repeat_filter(observations):
  '''
  Takes Python array of Pfam IDs in ant filters it for repeats.
  '''
  global repeats

  return [o for o in observations if o not in repeats]
