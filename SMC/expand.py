#!/usr/bin/env python

import os
import shutil
import sys
import getopt
import string
import re

def expand(forg):
    bas, ext = os.path.splitext(forg)
    fopt = bas+'_exp'+ext

    f0 = open(forg,'r')
    f1 = open(fopt,'w')

    contRE = re.compile(r"&\s*\n")

    multiline = []
    for line in f0.readlines():

        if fortran_comment(line):
            f1.write(line)
            continue 

        if contRE.search(line):

            multiline.append(line)

        elif multiline:  # the last line of multiline

            multiline.append(line)

            singleline = merge_lines(multiline)
            newMultiLine = expand_line(singleline)
            if newMultiLine:
                for l in multiline:
                    f1.write(comment_out(l))                
                f1.write(newMultiLine)
            else:
                for l in multiline:
                    f1.write(l)

            multiline = []

        else:

            newline = expand_line(line)
            if newline:
                f1.write(comment_out(line))
                f1.write(newline)
            else:
                f1.write(line)

    f1.close()
    f0.close()


def expand_line(line):
    mmRE = re.compile(r"matmul\s*\(.*\)")
    dpRE = re.compile(r"dot_product\s*\(.*matmul\s*\(.*\).*\)")
    search_mm = mmRE.search(line)
    search_dp = dpRE.search(line)
    if search_mm and not search_dp:
        return expand_matmul(line)
    elif search_dp:
        return expand_aMu(line)
    else:
        return ''

def expand_matmul(line):
    m8RE = re.compile(r"matmul\s*\(.*M8.*\)")
    if m8RE.search(line):
        return expand_matmul_8(line)
    else:
        return ''

def expand_aMu(line):
    m8RE = re.compile(r"matmul\s*\(.*M8.*\)")
    if m8RE.search(line):
        return expand_aMu_8(line)
    else:
        return ''

def expand_matmul_8(line):
    # expand (1) lhs = matmul(M8, u( .... )) OR
    #        (2) lhs = matmul(a( .... ), M8)
    lhs, rhs = line.split('=')

    lhs = lhs.strip(' \t')
    i = string.find(line, lhs)
    indent = line[0:i]
    moreindent = indent+'   '
    for i in range(len(lhs)):
        moreindent = moreindent+' '

    rhs = rhs.strip(' \t\n\r')
    args = rhs[7:-1].replace(' ','')  # 7 comes form 'matmul('
    if args[0:2] == 'M8':  
        # M8,u(,,,)
        x = args[3:]
        u = expand_fortran_slice(x)
        return indent+lhs+'(1) = '+'M8(1,1) * '+u[0]+' &\n' + \
               moreindent + ' + ' +'M8(1,2) * '+u[1]+' &\n' + \
               moreindent + ' + ' +'M8(1,3) * '+u[2]+' &\n' + \
               moreindent + ' + ' +'M8(1,4) * '+u[3]+' &\n' + \
               moreindent + ' - ' +'M8(8,4) * '+u[4]+'\n'   + \
               indent+lhs+'(2) = '+'M8(2,1) * '+u[0]+' &\n' + \
               moreindent + ' + ' +'M8(2,2) * '+u[1]+' &\n' + \
               moreindent + ' + ' +'M8(2,3) * '+u[2]+' &\n' + \
               moreindent + ' + ' +'M8(2,4) * '+u[3]+' &\n' + \
               moreindent + ' - ' +'M8(7,4) * '+u[4]+' &\n' + \
               moreindent + ' - ' +'M8(7,3) * '+u[5]+'\n'   + \
               indent+lhs+'(3) = '+'M8(3,1) * '+u[0]+' &\n' + \
               moreindent + ' + ' +'M8(3,2) * '+u[1]+' &\n' + \
               moreindent + ' + ' +'M8(3,3) * '+u[2]+' &\n' + \
               moreindent + ' + ' +'M8(3,4) * '+u[3]+' &\n' + \
               moreindent + ' - ' +'M8(6,4) * '+u[4]+' &\n' + \
               moreindent + ' - ' +'M8(6,3) * '+u[5]+' &\n' + \
               moreindent + ' - ' +'M8(6,2) * '+u[6]+'\n'   + \
               indent+lhs+'(4) = '+'M8(4,1) * '+u[0]+' &\n' + \
               moreindent + ' + ' +'M8(4,2) * '+u[1]+' &\n' + \
               moreindent + ' + ' +'M8(4,3) * '+u[2]+' &\n' + \
               moreindent + ' + ' +'M8(4,4) * '+u[3]+' &\n' + \
               moreindent + ' - ' +'M8(5,4) * '+u[4]+' &\n' + \
               moreindent + ' - ' +'M8(5,3) * '+u[5]+' &\n' + \
               moreindent + ' - ' +'M8(5,2) * '+u[6]+' &\n' + \
               moreindent + ' - ' +'M8(5,1) * '+u[7]+'\n'   + \
               indent+lhs+'(5) = '+'M8(5,1) * '+u[0]+' &\n' + \
               moreindent + ' + ' +'M8(5,2) * '+u[1]+' &\n' + \
               moreindent + ' + ' +'M8(5,3) * '+u[2]+' &\n' + \
               moreindent + ' + ' +'M8(5,4) * '+u[3]+' &\n' + \
               moreindent + ' - ' +'M8(4,4) * '+u[4]+' &\n' + \
               moreindent + ' - ' +'M8(4,3) * '+u[5]+' &\n' + \
               moreindent + ' - ' +'M8(4,2) * '+u[6]+' &\n' + \
               moreindent + ' - ' +'M8(4,1) * '+u[7]+'\n'   + \
               indent+lhs+'(6) = '+'M8(6,2) * '+u[1]+' &\n' + \
               moreindent + ' + ' +'M8(6,3) * '+u[2]+' &\n' + \
               moreindent + ' + ' +'M8(6,4) * '+u[3]+' &\n' + \
               moreindent + ' - ' +'M8(3,4) * '+u[4]+' &\n' + \
               moreindent + ' - ' +'M8(3,3) * '+u[5]+' &\n' + \
               moreindent + ' - ' +'M8(3,2) * '+u[6]+' &\n' + \
               moreindent + ' - ' +'M8(3,1) * '+u[7]+'\n'   + \
               indent+lhs+'(7) = '+'M8(7,3) * '+u[2]+' &\n' + \
               moreindent + ' + ' +'M8(7,4) * '+u[3]+' &\n' + \
               moreindent + ' - ' +'M8(2,4) * '+u[4]+' &\n' + \
               moreindent + ' - ' +'M8(2,3) * '+u[5]+' &\n' + \
               moreindent + ' - ' +'M8(2,2) * '+u[6]+' &\n' + \
               moreindent + ' - ' +'M8(2,1) * '+u[7]+'\n'   + \
               indent+lhs+'(8) = '+'M8(8,4) * '+u[3]+' &\n' + \
               moreindent + ' - ' +'M8(1,4) * '+u[4]+' &\n' + \
               moreindent + ' - ' +'M8(1,3) * '+u[5]+' &\n' + \
               moreindent + ' - ' +'M8(1,2) * '+u[6]+' &\n' + \
               moreindent + ' - ' +'M8(1,1) * '+u[7]+'\n' 
    else: 
        # a(,,,), M8
        x = args[0:-3]
        a = expand_fortran_slice(x)
        return indent+lhs+'(1) = '+a[0]+' * M8(1,1) &\n' + \
               moreindent + ' + ' +a[1]+' * M8(2,1) &\n' + \
               moreindent + ' + ' +a[2]+' * M8(3,1) &\n' + \
               moreindent + ' + ' +a[3]+' * M8(4,1) &\n' + \
               moreindent + ' + ' +a[4]+' * M8(5,1)\n'   + \
               indent+lhs+'(2) = '+a[0]+' * M8(1,2) &\n' + \
               moreindent + ' + ' +a[1]+' * M8(2,2) &\n' + \
               moreindent + ' + ' +a[2]+' * M8(3,2) &\n' + \
               moreindent + ' + ' +a[3]+' * M8(4,2) &\n' + \
               moreindent + ' + ' +a[4]+' * M8(5,2) &\n' + \
               moreindent + ' + ' +a[5]+' * M8(6,2)\n'   + \
               indent+lhs+'(3) = '+a[0]+' * M8(1,3) &\n' + \
               moreindent + ' + ' +a[1]+' * M8(2,3) &\n' + \
               moreindent + ' + ' +a[2]+' * M8(3,3) &\n' + \
               moreindent + ' + ' +a[3]+' * M8(4,3) &\n' + \
               moreindent + ' + ' +a[4]+' * M8(5,3) &\n' + \
               moreindent + ' + ' +a[5]+' * M8(6,3) &\n' + \
               moreindent + ' + ' +a[6]+' * M8(7,3)\n'   + \
               indent+lhs+'(4) = '+a[0]+' * M8(1,4) &\n' + \
               moreindent + ' + ' +a[1]+' * M8(2,4) &\n' + \
               moreindent + ' + ' +a[2]+' * M8(3,4) &\n' + \
               moreindent + ' + ' +a[3]+' * M8(4,4) &\n' + \
               moreindent + ' + ' +a[4]+' * M8(5,4) &\n' + \
               moreindent + ' + ' +a[5]+' * M8(6,4) &\n' + \
               moreindent + ' + ' +a[6]+' * M8(7,4) &\n' + \
               moreindent + ' + ' +a[7]+' * M8(8,4)\n'   + \
               indent+lhs+'(5) =-'+a[0]+' * M8(8,4) &\n' + \
               moreindent + ' - ' +a[1]+' * M8(7,4) &\n' + \
               moreindent + ' - ' +a[2]+' * M8(6,4) &\n' + \
               moreindent + ' - ' +a[3]+' * M8(5,4) &\n' + \
               moreindent + ' - ' +a[4]+' * M8(4,4) &\n' + \
               moreindent + ' - ' +a[5]+' * M8(3,4) &\n' + \
               moreindent + ' - ' +a[6]+' * M8(2,4) &\n' + \
               moreindent + ' - ' +a[7]+' * M8(1,4)\n'   + \
               indent+lhs+'(6) =-'+a[1]+' * M8(7,3) &\n' + \
               moreindent + ' - ' +a[2]+' * M8(6,3) &\n' + \
               moreindent + ' - ' +a[3]+' * M8(5,3) &\n' + \
               moreindent + ' - ' +a[4]+' * M8(4,3) &\n' + \
               moreindent + ' - ' +a[5]+' * M8(3,3) &\n' + \
               moreindent + ' - ' +a[6]+' * M8(2,3) &\n' + \
               moreindent + ' - ' +a[7]+' * M8(1,3)\n' + \
               indent+lhs+'(7) =-'+a[2]+' * M8(6,2) &\n' + \
               moreindent + ' - ' +a[3]+' * M8(5,2) &\n' + \
               moreindent + ' - ' +a[4]+' * M8(4,2) &\n' + \
               moreindent + ' - ' +a[5]+' * M8(3,2) &\n' + \
               moreindent + ' - ' +a[6]+' * M8(2,2) &\n' + \
               moreindent + ' - ' +a[7]+' * M8(1,2)\n' + \
               indent+lhs+'(8) =-'+a[3]+' * M8(5,1) &\n' + \
               moreindent + ' - ' +a[4]+' * M8(4,1) &\n' + \
               moreindent + ' - ' +a[5]+' * M8(3,1) &\n' + \
               moreindent + ' - ' +a[6]+' * M8(2,1) &\n' + \
               moreindent + ' - ' +a[7]+' * M8(1,1)\n' 


def expand_aMu_8(line):
    # expand lhs = dot_product(matmul(a,M8),u) [+ dot_product()]
    lhs, rhs = line.split('=')

    lhs = lhs.strip(' \t')
    i = string.find(line, lhs)
    indent = line[0:i]
    moreindent = indent+'   '

    rhs = rhs.strip(' \t\n\r').replace(' ','')
    rhs = rhs[19:] # strip dot_product(matmul(
    i = string.find(rhs,')+dot_product')
    if i > 0:
        aMu = rhs[0:i]
        dp2 = rhs[i+2:]
        line0 = indent+lhs+' = '+dp2+' + &\n'
    else:
        aMu = rhs[:-1]
        line0 = indent+lhs+' = &\n'

    a, u = aMu.split(',M8),')
    aL = expand_fortran_slice(a)
    uL = expand_fortran_slice(u)

    return line0 + \
           moreindent + '( M8(1,1)*('+aL[0]+'*'+uL[0]+'-'+aL[7]+'*'+uL[7]+') &\n' +\
           moreindent + '+ M8(2,1)*('+aL[1]+'*'+uL[0]+'-'+aL[6]+'*'+uL[7]+') &\n' +\
           moreindent + '+ M8(3,1)*('+aL[2]+'*'+uL[0]+'-'+aL[5]+'*'+uL[7]+') &\n' +\
           moreindent + '+ M8(4,1)*('+aL[3]+'*'+uL[0]+'-'+aL[4]+'*'+uL[7]+') &\n' +\
           moreindent + '+ M8(5,1)*('+aL[4]+'*'+uL[0]+'-'+aL[3]+'*'+uL[7]+') &\n' +\
           \
           moreindent + '+ M8(1,2)*('+aL[0]+'*'+uL[1]+'-'+aL[7]+'*'+uL[6]+') &\n' +\
           moreindent + '+ M8(2,2)*('+aL[1]+'*'+uL[1]+'-'+aL[6]+'*'+uL[6]+') &\n' +\
           moreindent + '+ M8(3,2)*('+aL[2]+'*'+uL[1]+'-'+aL[5]+'*'+uL[6]+') &\n' +\
           moreindent + '+ M8(4,2)*('+aL[3]+'*'+uL[1]+'-'+aL[4]+'*'+uL[6]+') &\n' +\
           moreindent + '+ M8(5,2)*('+aL[4]+'*'+uL[1]+'-'+aL[3]+'*'+uL[6]+') &\n' +\
           moreindent + '+ M8(6,2)*('+aL[5]+'*'+uL[1]+'-'+aL[2]+'*'+uL[6]+') &\n' +\
           \
           moreindent + '+ M8(1,3)*('+aL[0]+'*'+uL[2]+'-'+aL[7]+'*'+uL[5]+') &\n' +\
           moreindent + '+ M8(2,3)*('+aL[1]+'*'+uL[2]+'-'+aL[6]+'*'+uL[5]+') &\n' +\
           moreindent + '+ M8(3,3)*('+aL[2]+'*'+uL[2]+'-'+aL[5]+'*'+uL[5]+') &\n' +\
           moreindent + '+ M8(4,3)*('+aL[3]+'*'+uL[2]+'-'+aL[4]+'*'+uL[5]+') &\n' +\
           moreindent + '+ M8(5,3)*('+aL[4]+'*'+uL[2]+'-'+aL[3]+'*'+uL[5]+') &\n' +\
           moreindent + '+ M8(6,3)*('+aL[5]+'*'+uL[2]+'-'+aL[2]+'*'+uL[5]+') &\n' +\
           moreindent + '+ M8(7,3)*('+aL[6]+'*'+uL[2]+'-'+aL[1]+'*'+uL[5]+') &\n' +\
           \
           moreindent + '+ M8(1,4)*('+aL[0]+'*'+uL[3]+'-'+aL[7]+'*'+uL[4]+') &\n' +\
           moreindent + '+ M8(2,4)*('+aL[1]+'*'+uL[3]+'-'+aL[6]+'*'+uL[4]+') &\n' +\
           moreindent + '+ M8(3,4)*('+aL[2]+'*'+uL[3]+'-'+aL[5]+'*'+uL[4]+') &\n' +\
           moreindent + '+ M8(4,4)*('+aL[3]+'*'+uL[3]+'-'+aL[4]+'*'+uL[4]+') &\n' +\
           moreindent + '+ M8(5,4)*('+aL[4]+'*'+uL[3]+'-'+aL[3]+'*'+uL[4]+') &\n' +\
           moreindent + '+ M8(6,4)*('+aL[5]+'*'+uL[3]+'-'+aL[2]+'*'+uL[4]+') &\n' +\
           moreindent + '+ M8(7,4)*('+aL[6]+'*'+uL[3]+'-'+aL[1]+'*'+uL[4]+') &\n' +\
           moreindent + '+ M8(8,4)*('+aL[7]+'*'+uL[3]+'-'+aL[0]+'*'+uL[4]+') )\n' 


def expand_fortran_slice(x):
    xs = []
    lb = string.find(x, '(')
    rb = string.find(x, ')')
    v = x[0:lb]
    indices = x[lb+1:rb].split(',')
    i_index = indices[0]
    j_index = indices[1]
    k_index = indices[2]
    if string.find(i_index,':') >= 0:
        istart = int(i_index.split(':')[0][1:])
        for i in range(istart,-istart):
            xs.append(v+'(i+'+str(i)+',j,k,NCOMP)')
    elif string.find(j_index,':') >= 0:
        jstart = int(j_index.split(':')[0][1:])
        for j in range(jstart,-jstart):
            xs.append(v+'(i,j+'+str(j)+',k,NCOMP)')
    else:
        kstart = int(k_index.split(':')[0][1:])
        for k in range(kstart,-kstart):
            xs.append(v+'(i,j,k+'+str(k)+',NCOMP)')
    if len(indices) == 4:
        return [t.replace('NCOMP',indices[3]).replace('+-','-').replace('+0','  ') for t in xs]
    else:
        return [t.replace(',NCOMP','').replace('+-','-').replace('+0','  ') for t in xs]

def merge_lines(multiline):
    lns = []
    lns.append(multiline[0].rstrip(' \t\n\r&'))
    for l in multiline[1:]:
        lns.append(l.strip(' \t\n\r&'))
    return ' '.join(lns)+'\n'

def comment_out(line):
    return '!'+line

def fortran_comment(line):
    commRE = re.compile(r"\s*!")
    if commRE.match(line):
        return True
    else:
        return False

if __name__== "__main__":
    if len(sys.argv) == 1:
        print "usage: expand.py kernels.f90"
        sys.exit(1)
    expand(sys.argv[1])
