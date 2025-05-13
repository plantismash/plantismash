# vim: set fileencoding=utf-8 :
#
# Copyright (C) 2010-2012 Marnix H. Medema
# University of Groningen
# Department of Microbial Physiology / Groningen Bioinformatics Centre
#
# Copyright (C) 2011,2012 Kai Blin
# University of Tuebingen
# Interfaculty Institute of Microbiology and Infection Medicine
# Div. of Microbiology/Biotechnology
#
# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

import os
import logging
import random
from antismash import utils
from antismash.lib.pysvg.filter import *
from antismash.lib.pysvg.gradient import *
from antismash.lib.pysvg.linking import *
from antismash.lib.pysvg.script import *
from antismash.lib.pysvg.shape import *
from antismash.lib.pysvg.structure import *
from antismash.lib.pysvg.style import *
from antismash.lib.pysvg.text import *
from antismash.lib.pysvg.builders import *
import pprint

def create_svgs(options, seq_record):
    logging.info("Writing visualization SVGs")
    screenwidth = 1024
    seq_record.geneposdict = create_genecluster_svg(options, seq_record, screenwidth)
    if options.clusterblast:
        seq_record.clusterblastpositiondata = create_clusterblast_svg(options, seq_record, screenwidth)
    if options.subclusterblast:
        seq_record.subclusterblastpositiondata = create_clusterblast_svg(options, seq_record, screenwidth, searchtype="subclusters")
    if options.knownclusterblast:
        seq_record.knownclusterblastpositiondata = create_clusterblast_svg(options, seq_record, screenwidth, searchtype="knownclusters")

def create_genecluster_svg(options, seq_record, screenwidth):
    #Create genecluster svg for each gene cluster
    geneposdict = {}
    geneclusters = utils.get_sorted_cluster_features(seq_record)
    for genecluster in geneclusters:
        qclusternr = utils.get_cluster_number(genecluster)
        data = seq_record.qgeneclusterdata[qclusternr]
        clustersize = data[1]
        genes = data[2]
        starts = data[4]
        ends = data[5]
        strands = data[6]
        pksnrpsprots = data[7]
        pksnrpsdomains = data[9]
        colors = data[16]
        relpositions = relativepositions(starts,ends,clustersize, screenwidth)
        rel_starts = relpositions[0]
        rel_ends = relpositions[1]
        y = 0
        for i in genes:
            if int(starts[y]) < int(ends[y]):
                geneposdict[i] = [int(starts[y]),int(ends[y])]
            else:
                geneposdict[i] = [int(ends[y]),int(starts[y])]
            y += 1
        s = geneclustersvg(genes, rel_starts,rel_ends, strands, geneposdict, pksnrpsprots, pksnrpsdomains, qclusternr, colors, screenwidth)
        outfile = open(options.svgdir + os.sep + "genecluster" + str(qclusternr) + ".svg","w")
        outfile.write(s.getXML())
        outfile.close()
    return geneposdict

def _gene_arrow(start,end,strand,color,base,height):
    halfheight = height/2
    if start > end:
        start2 = end
        end2 = start
        start = start2
        end = end2
    oh = ShapeBuilder()
    if (end - start) < halfheight:
        if (strand == "+"):
            pointsAsTuples=[(start,base),
                            (end,base - halfheight),
                            (start,base - height),
                            (start,base)
                           ]
        if (strand == "-"):
            pointsAsTuples=[(start,base - halfheight),
                            (end,base - height),
                            (end,base),
                            (start,base - halfheight)
                           ]
    else:
        if (strand == "+"):
            arrowstart = end-halfheight
            pointsAsTuples=[(start,base),
                            (arrowstart,base),
                            (end,base-halfheight),
                            (arrowstart,base - height),
                            (start,base - height),
                            (start,base)
                           ]
        if (strand == "-"):
            arrowstart = start + halfheight
            pointsAsTuples=[(start,base - halfheight),
                            (arrowstart,base - height),
                            (end,base - height),
                            (end,base),
                            (arrowstart,base),
                            (start,base - halfheight)
                           ]
    pg=oh.createPolygon(points=oh.convertTupleArrayToPoints(pointsAsTuples),strokewidth=1, stroke='black', fill=color)
    return pg

def generate_rgbscheme(nr):
    usablenumbers = [1,2,4,8,12,18,24,32,48,64,10000]
    lengthsdict = {1:[1,1,1],2:[1,1,2],4:[1,2,2],8:[2,2,2],12:[2,2,3],18:[2,3,3],24:[3,3,3],32:[3,3,4],48:[3,4,4],64:[4,4,4]}
    shortestdistance = 10000
    for i in usablenumbers:
        distance = i - nr
        if distance >= 0:
            if distance < shortestdistance:
                shortestdistance = distance
                closestnr = i
    toohigh = "n"
    if closestnr == 10000:
        toohigh = "y"
        closestnr = 64
    xyznumbers = lengthsdict[closestnr]
    x = xyznumbers[0]
    y = xyznumbers[1]
    z = xyznumbers[2]
    xpoints = []
    xpoint = (255/z)/2
    for i in range(x):
        xpoints.append(xpoint)
        xpoint += (255/x)
    ypoints = []
    ypoint = (255/z)/2
    for i in range(y):
        ypoints.append(ypoint)
        ypoint += (255/y)
    zpoints = []
    zpoint = (255/z)/2
    for i in range(z):
        zpoints.append(zpoint)
        zpoint += (255/z)
    colorlist = []
    for i in xpoints:
        for j in ypoints:
            for k in zpoints:
                rgb = "rgb(" + str(i) + "," + str(j) + "," + str(k) + ")"
                colorlist.append(rgb)
    if toohigh == "y":
        colorlist = colorlist + colorlist + colorlist + colorlist + colorlist + colorlist + colorlist + colorlist + colorlist + colorlist + colorlist + colorlist + colorlist + colorlist + colorlist + colorlist + colorlist + colorlist + colorlist + colorlist
    if closestnr == 24:
        colorlist = colorlist[:15] + colorlist[18:]
    if closestnr == 32:
        colorlist = colorlist[:21] + colorlist[24:]
    colorlist2 = []
    if closestnr == 1:
        colorlist2.append("red")
    if closestnr == 2:
        colorlist2.append("red")
        colorlist2.append("green")
    if closestnr == 4:
        colorlist2.append("red")
        colorlist2.append("green")
        colorlist2.append("blue")
        colorlist2.append("yellow")
    if closestnr == 8:
        neworder=[4,1,2,5,6,7,3,0]
        colorlist2 = [colorlist[i] for i in neworder]
    if closestnr == 12:
        neworder=[6,3,5,9,7,2,11,4,8,1,10,0]
        colorlist2 = [colorlist[i] for i in neworder]
    if closestnr == 18:
        neworder=[9,6,2,14,15,8,12,10,3,5,7,11,4,1,16,13,0]
        colorlist2 = [colorlist[i] for i in neworder]
    if closestnr == 24:
        neworder=[15,12,9,6,5,0,21,1,16,14,8,17,2,23,22,3,13,7,10,4,18,20,19,11]
        colorlist2 = [colorlist[i] for i in neworder]
    if closestnr == 32:
        neworder = [21,19,27,6,8,1,14,7,20,13,9,30,4,23,18,12,5,29,24,17,11,31,2,28,22,15,26,3,20,16,10,25]
        colorlist2 = [colorlist[i] for i in neworder]
    if closestnr > 32:
        random.shuffle(colorlist)
        colorlist2 = colorlist
    colorlist = colorlist2
    return colorlist

def startendsitescheck(starts,ends):
    #Check whether start sites are always lower than end sites, reverse if necessary
    starts2 = []
    ends2 = []
    a = 0
    for i in starts:
        if int(i) > int(ends[a]):
            starts2.append(ends[a])
            ends2.append(i)
        else:
            starts2.append(i)
            ends2.append(ends[a])
        a += 1
    ends = ends2
    starts = starts2
    return [starts,ends]

def LinearGradient(startcolor,stopcolor,gradientname):
    d = defs()
    lg = linearGradient()
    lg.set_id(gradientname)
    s = stop(offset="0%")
    s.set_stop_color(startcolor)
    s.set_stop_opacity(1)
    lg.addElement(s)
    s = stop(offset="100%")
    s.set_stop_color(stopcolor)
    s.set_stop_opacity(1)
    lg.addElement(s)
    d.addElement(lg)
    return d

def domains_svg_part(s, geneposdict, pksnrpsprots, pksnrpsdomains, qclusternr, screenwidth):
    #Add domain depictions to SVG
    oh = ShapeBuilder()
    group = g()
    #Determine longest protein to decide on scaling
    longestprot = 1
    protlengthdict = {}
    for i in pksnrpsprots:
        protlength = (max(geneposdict[i]) - min(geneposdict[i])) / 3
        protlengthdict[i] = protlength
        if protlength > longestprot:
            longestprot = protlength
    z = 0
    w = 0
    ksnr = 1
    atnr = 1
    dhnr = 1
    krnr = 1
    ernr = 1
    acpnr = 1
    cnr = 1
    enr = 1
    anr = 1
    pcpnr = 1
    tenr = 1
    othernr = 1
    for i in pksnrpsprots:
        domains = pksnrpsdomains[i][0]
        domainsdict = pksnrpsdomains[i][1]
        protlength = protlengthdict[i]
        group.addElement(oh.createLine(10,(z * 60 ),10 + ((float(protlength) / float(longestprot)) * (screenwidth * 0.75)),(z * 60 ), strokewidth = 1, stroke = "grey"))
        s.addElement(group)
        aa2pixelratio = float(longestprot / float(screenwidth * 0.75))
        myStyle = StyleBuilder()
        myStyle.setFontFamily(fontfamily="MS Reference Sans Serif")
        myStyle.setFontWeight(fontweight='bold')
        myStyle.setFontSize('12px')
        for j in domains:
            startpos = domainsdict[j][0]
            endpos = domainsdict[j][1]
            if "PKS_KS" in j:
                c = LinearGradient("#08B208","#81F781","KS_domain" + str(qclusternr) + "_" + str(ksnr))
                d = LinearGradient("#81F781","#08B208","KS_line" + str(qclusternr) + "_" + str(ksnr))
                e = oh.createRect(str(10 + startpos / aa2pixelratio),str((z * 60 ) - 8),str((endpos-startpos) / aa2pixelratio),15,8,strokewidth=1,stroke='url(#KS_line' + str(qclusternr) + "_" + str(ksnr) + ")",fill="url(#KS_domain" + str(qclusternr) + "_" + str(ksnr) + ")")
                f = text("KS",((-4 + startpos / aa2pixelratio) + 0.5 * ((endpos-startpos) / aa2pixelratio)), ((z * 60 ) + 4),fill='#0A2A0A')
                if ((endpos-startpos) / aa2pixelratio) < 100 and ((endpos-startpos) / aa2pixelratio) >= 20:
                    myStyle.setFontSize('8px')
                    f = text("KS",((startpos / aa2pixelratio) + 0.5 * ((endpos-startpos) / aa2pixelratio)), ((z * 60 ) + 3),fill='#3B0B0B')
                elif ((endpos-startpos) / aa2pixelratio) < 20:
                    f = "notext"
                if f != "notext":
                    f.set_style(myStyle.getStyle())
                myStyle.setFontSize('12px')
                group = g()
                group.addElement(c)
                group.addElement(d)
                group.addElement(e)
                if f != "notext":
                    group.addElement(f)
                group.set_id("b" + str(qclusternr) + "_00%s"%w)
                s.addElement(group)
                ksnr += 1
            elif "PKS_AT" in j:
                c = LinearGradient("#DC0404","#F78181","AT_domain" + str(qclusternr) + "_" + str(atnr))
                d = LinearGradient("#F78181","#DC0404","AT_line" + str(qclusternr) + "_" + str(atnr))
                e = oh.createRect(str(10 + startpos / aa2pixelratio),str((z * 60 ) - 8),str((endpos-startpos) / aa2pixelratio),15,8,strokewidth=1,stroke='url(#AT_line' + str(qclusternr) + "_" + str(atnr) + ")",fill="url(#AT_domain" + str(qclusternr) + "_" + str(atnr) + ")")
                f = text("AT",((-4 + startpos / aa2pixelratio) + 0.5 * ((endpos-startpos) / aa2pixelratio)), ((z * 60 ) + 4),fill='#2A1B0A')
                if ((endpos-startpos) / aa2pixelratio) < 100 and ((endpos-startpos) / aa2pixelratio) >= 20:
                    myStyle.setFontSize('8px')
                    f = text("AT",((startpos / aa2pixelratio) + 0.5 * ((endpos-startpos) / aa2pixelratio)), ((z * 60 ) + 3),fill='#2A1B0A')
                elif ((endpos-startpos) / aa2pixelratio) < 20:
                    f = "notext"
                if f != "notext":
                    f.set_style(myStyle.getStyle())
                myStyle.setFontSize('12px')
                group = g()
                group.addElement(c)
                group.addElement(d)
                group.addElement(e)
                if f != "notext":
                    group.addElement(f)
                group.set_id("b" + str(qclusternr) + "_00%s"%w)
                s.addElement(group)
                atnr += 1
            elif "PKS_DH" in j or "PKS_DH2" in j:
                c = LinearGradient("#B45F04","#F7BE81","DH_domain" + str(qclusternr) + "_" + str(dhnr))
                d = LinearGradient("#F7BE81","#B45F04","DH_line" + str(qclusternr) + "_" + str(dhnr))
                e = oh.createRect(str(10 + startpos / aa2pixelratio),str((z * 60 ) - 8),str((endpos-startpos) / aa2pixelratio),15,8,strokewidth=1,stroke='url(#DH_line' + str(qclusternr) + "_" + str(dhnr) + ")",fill="url(#DH_domain" + str(qclusternr) + "_" + str(dhnr) + ")")
                f = text("DH",((-4 + startpos / aa2pixelratio) + 0.5 * ((endpos-startpos) / aa2pixelratio)), ((z * 60 ) + 4),fill='#3B0B0B')
                if ((endpos-startpos) / aa2pixelratio) < 100 and ((endpos-startpos) / aa2pixelratio) >= 20:
                    myStyle.setFontSize('8px')
                    f = text("DH",((startpos / aa2pixelratio) + 0.5 * ((endpos-startpos) / aa2pixelratio)), ((z * 60 ) + 3),fill='#3B0B0B')
                elif ((endpos-startpos) / aa2pixelratio) < 20:
                    f = "notext"
                if f != "notext":
                    f.set_style(myStyle.getStyle())
                myStyle.setFontSize('12px')
                group = g()
                group.addElement(c)
                group.addElement(d)
                group.addElement(e)
                if f != "notext":
                    group.addElement(f)
                group.set_id("b" + str(qclusternr) + "_00%s"%w)
                s.addElement(group)
                dhnr += 1
            elif "PKS_KR" in j:
                c = LinearGradient("#089E4B","#81F781","KR_domain" + str(qclusternr) + "_" + str(krnr))
                d = LinearGradient("#81F781","#089E4B","KR_line" + str(qclusternr) + "_" + str(krnr))
                e = oh.createRect(str(10 + startpos / aa2pixelratio),str((z * 60 ) - 8),str((endpos-startpos) / aa2pixelratio),15,8,strokewidth=1,stroke='url(#KR_line' + str(qclusternr) + "_" + str(krnr) + ")",fill="url(#KR_domain" + str(qclusternr) + "_" + str(krnr) + ")")
                f = text("KR",((-4 + startpos / aa2pixelratio) + 0.5 * ((endpos-startpos) / aa2pixelratio)), ((z * 60 ) + 4),fill='#0A2A1B')
                if ((endpos-startpos) / aa2pixelratio) < 100 and ((endpos-startpos) / aa2pixelratio) >= 20:
                    myStyle.setFontSize('8px')
                    f = text("KR",((startpos / aa2pixelratio) + 0.5 * ((endpos-startpos) / aa2pixelratio)), ((z * 60 ) + 3),fill='#0A2A1B')
                elif ((endpos-startpos) / aa2pixelratio) < 20:
                    f = "notext"
                if f != "notext":
                    f.set_style(myStyle.getStyle())
                myStyle.setFontSize('12px')
                group = g()
                group.addElement(c)
                group.addElement(d)
                group.addElement(e)
                if f != "notext":
                    group.addElement(f)
                group.set_id("b" + str(qclusternr) + "_00%s"%w)
                s.addElement(group)
                krnr += 1
            elif "PKS_ER" in j:
                c = LinearGradient("#089E85","#81F7F3","ER_domain" + str(qclusternr) + "_" + str(ernr))
                d = LinearGradient("#81F7F3","#089E85","ER_line" + str(qclusternr) + "_" + str(ernr))
                e = oh.createRect(str(10 + startpos / aa2pixelratio),str((z * 60 ) - 8),str((endpos-startpos) / aa2pixelratio),15,8,strokewidth=1,stroke='url(#ER_line' + str(qclusternr) + "_" + str(ernr) + ")",fill="url(#ER_domain" + str(qclusternr) + "_" + str(ernr) + ")")
                f = text("ER",((-4 + startpos / aa2pixelratio) + 0.5 * ((endpos-startpos) / aa2pixelratio)), ((z * 60 ) + 4),fill='#0A2A29')
                if ((endpos-startpos) / aa2pixelratio) < 100 and ((endpos-startpos) / aa2pixelratio) >= 20:
                    myStyle.setFontSize('8px')
                    f = text("ER",((startpos / aa2pixelratio) + 0.5 * ((endpos-startpos) / aa2pixelratio)), ((z * 60 ) + 3),fill='#0A2A29')
                elif ((endpos-startpos) / aa2pixelratio) < 20:
                    f = "notext"
                if f != "notext":
                    f.set_style(myStyle.getStyle())
                myStyle.setFontSize('12px')
                group = g()
                group.addElement(c)
                group.addElement(d)
                group.addElement(e)
                if f != "notext":
                    group.addElement(f)
                group.set_id("b" + str(qclusternr) + "_00%s"%w)
                s.addElement(group)
                ernr += 1
            elif "ACP" in j:
                c = LinearGradient("#084BC6","#81BEF7","ACP_domain" + str(qclusternr) + "_" + str(acpnr))
                d = LinearGradient("#81BEF7","#084BC6","ACP_line" + str(qclusternr) + "_" + str(acpnr))
                e = oh.createRect(str(10 + startpos / aa2pixelratio),str((z * 60 ) - 8),str((endpos-startpos) / aa2pixelratio),15,8,strokewidth=1,stroke='url(#ACP_line' + str(qclusternr) + "_" + str(acpnr) + ")",fill="url(#ACP_domain" + str(qclusternr) + "_" + str(acpnr) + ")")
                f = text("ACP",((startpos / aa2pixelratio) + 0.5 * ((endpos-startpos) / aa2pixelratio)), ((z * 60 ) + 4),fill='#0A1B2A')
                if ((endpos-startpos) / aa2pixelratio) < 100 and ((endpos-startpos) / aa2pixelratio) >= 20:
                    myStyle.setFontSize('8px')
                    f = text("ACP",((-2 + startpos / aa2pixelratio) + 0.5 * ((endpos-startpos) / aa2pixelratio)), ((z * 60 ) + 3),fill='#0A1B2A')
                elif ((endpos-startpos) / aa2pixelratio) < 20:
                    f = "notext"
                if f != "notext":
                    f.set_style(myStyle.getStyle())
                myStyle.setFontSize('12px')
                group = g()
                group.addElement(c)
                group.addElement(d)
                group.addElement(e)
                if f != "notext":
                    group.addElement(f)
                group.set_id("b" + str(qclusternr) + "_00%s"%w)
                s.addElement(group)
                acpnr += 1
            elif ("C" in j or "Heterocyclization" in j ) and "ACP" not in j and "PCP" not in j and "NRPS-COM" not in j and "CAL" not in j and "cMT" not in j:
                c = LinearGradient("#393989","#8181F7","C_domain" + str(qclusternr) + "_" + str(cnr))
                d = LinearGradient("#8181F7","#393989","C_line" + str(qclusternr) + "_" + str(cnr))
                e = oh.createRect(str(10 + startpos / aa2pixelratio),str((z * 60 ) - 8),str((endpos-startpos) / aa2pixelratio),15,8,strokewidth=1,stroke='url(#C_line' + str(qclusternr) + "_" + str(cnr) + ")",fill="url(#C_domain" + str(qclusternr) + "_" + str(cnr) + ")")
                f = text("C",((startpos / aa2pixelratio) + 0.5 * ((endpos-startpos) / aa2pixelratio)), ((z * 60 ) + 4),fill='#0A0A2A')
                if ((endpos-startpos) / aa2pixelratio) < 100 and ((endpos-startpos) / aa2pixelratio) >= 20:
                    myStyle.setFontSize('8px')
                    f = text("C",((startpos / aa2pixelratio) + 0.5 * ((endpos-startpos) / aa2pixelratio)), ((z * 60 ) + 3),fill='#0A0A2A')
                elif ((endpos-startpos) / aa2pixelratio) < 20:
                    f = "notext"
                if f != "notext":
                    f.set_style(myStyle.getStyle())
                myStyle.setFontSize('12px')
                group = g()
                group.addElement(c)
                group.addElement(d)
                group.addElement(e)
                if f != "notext":
                    group.addElement(f)
                group.set_id("b" + str(qclusternr) + "_00%s"%w)
                s.addElement(group)
                cnr += 1
            elif "Epimerization" in j and "ER" not in j and "TE" not in j:
                c = LinearGradient("#393989","#8181F7","E_domain" + str(qclusternr) + "_" + str(enr))
                d = LinearGradient("#8181F7","#393989","E_line" + str(qclusternr) + "_" + str(enr))
                e = oh.createRect(str(10 + startpos / aa2pixelratio),str((z * 60 ) - 8),str((endpos-startpos) / aa2pixelratio),15,8,strokewidth=1,stroke='url(#E_line' + str(qclusternr) + "_" + str(enr) + ")",fill="url(#E_domain" + str(qclusternr) + "_" + str(enr) + ")")
                f = text("E",((startpos / aa2pixelratio) + 0.5 * ((endpos-startpos) / aa2pixelratio)), ((z * 60 ) + 4),fill='#0A0A2A')
                if ((endpos-startpos) / aa2pixelratio) < 100 and ((endpos-startpos) / aa2pixelratio) >= 20:
                    myStyle.setFontSize('8px')
                    f = text("E",((startpos / aa2pixelratio) + 0.5 * ((endpos-startpos) / aa2pixelratio)), ((z * 60 ) + 3),fill='#0A0A2A')
                elif ((endpos-startpos) / aa2pixelratio) < 20:
                    f = "notext"
                if f != "notext":
                    f.set_style(myStyle.getStyle())
                myStyle.setFontSize('12px')
                group = g()
                group.addElement(c)
                group.addElement(d)
                group.addElement(e)
                if f != "notext":
                    group.addElement(f)
                group.set_id("b" + str(qclusternr) + "_00%s"%w)
                s.addElement(group)
                enr += 1
            elif ("AMP" in j or "A-OX" in j):
                c = LinearGradient("#56157F","#BE81F7","A_domain" + str(qclusternr) + "_" + str(anr))
                d = LinearGradient("#BE81F7","#56157F","A_line" + str(qclusternr) + "_" + str(anr))
                e = oh.createRect(str(10 + startpos / aa2pixelratio),str((z * 60 ) - 8),str((endpos-startpos) / aa2pixelratio),15,8,strokewidth=1,stroke='url(#A_line' + str(qclusternr) + "_" + str(anr) + ")",fill="url(#A_domain" + str(qclusternr) + "_" + str(anr) + ")")
                f = text("A",((startpos / aa2pixelratio) + 0.5 * ((endpos-startpos) / aa2pixelratio)), ((z * 60 ) + 4),fill='#1B0A2A')
                if ((endpos-startpos) / aa2pixelratio) < 100 and ((endpos-startpos) / aa2pixelratio) >= 20:
                    myStyle.setFontSize('8px')
                    f = text("A",((startpos / aa2pixelratio) + 0.5 * ((endpos-startpos) / aa2pixelratio)), ((z * 60 ) + 3),fill='#1B0A2A')
                elif ((endpos-startpos) / aa2pixelratio) < 20:
                    f = "notext"
                if f != "notext":
                    f.set_style(myStyle.getStyle())
                myStyle.setFontSize('12px')
                group = g()
                group.addElement(c)
                group.addElement(d)
                group.addElement(e)
                if f != "notext":
                    group.addElement(f)
                group.set_id("b" + str(qclusternr) + "_00%s"%w)
                s.addElement(group)
                anr += 1
            elif "PCP" in j:
                c = LinearGradient("#084BC6","#81BEF7","PCP_domain" + str(qclusternr) + "_" + str(pcpnr))
                d = LinearGradient("#81BEF7","#084BC6","PCP_line" + str(qclusternr) + "_" + str(pcpnr))
                e = oh.createRect(str(10 + startpos / aa2pixelratio),str((z * 60 ) - 8),str((endpos-startpos) / aa2pixelratio),15,8,strokewidth=1,stroke='url(#PCP_line' + str(qclusternr) + "_" + str(pcpnr) + ")",fill="url(#PCP_domain" + str(qclusternr) + "_" + str(pcpnr) + ")")
                f = text("PCP",((startpos / aa2pixelratio) + 0.5 * ((endpos-startpos) / aa2pixelratio)), ((z * 60 ) + 4),fill='#0A1B2A')
                if ((endpos-startpos) / aa2pixelratio) < 100 and ((endpos-startpos) / aa2pixelratio) >= 20:
                    myStyle.setFontSize('8px')
                    f = text("PCP",((-2 + startpos / aa2pixelratio) + 0.5 * ((endpos-startpos) / aa2pixelratio)), ((z * 60 ) + 3),fill='#0A1B2A')
                elif ((endpos-startpos) / aa2pixelratio) < 20:
                    f = "notext"
                if f != "notext":
                    f.set_style(myStyle.getStyle())
                myStyle.setFontSize('12px')
                group = g()
                group.addElement(c)
                group.addElement(d)
                group.addElement(e)
                if f != "notext":
                    group.addElement(f)
                group.set_id("b" + str(qclusternr) + "_00%s"%w)
                s.addElement(group)
                pcpnr += 1
            elif "Thioesterase" in j or "TD" in j:
                c = LinearGradient("#750072","#F5A9F2","TE_domain" + str(qclusternr) + "_" + str(tenr))
                d = LinearGradient("#F5A9F2","#750072","TE_line" + str(qclusternr) + "_" + str(tenr))
                e = oh.createRect(str(10 + startpos / aa2pixelratio),str((z * 60 ) - 8),str((endpos-startpos) / aa2pixelratio),15,8,strokewidth=1,stroke='url(#TE_line' + str(qclusternr) + "_" + str(tenr) + ")",fill="url(#TE_domain" + str(qclusternr) + "_" + str(tenr) + ")")
                if "Thioesterase" in j:
                    f = text("TE",((startpos / aa2pixelratio) + 0.5 * ((endpos-startpos) / aa2pixelratio)), ((z * 60 ) + 4),fill='#2A0A29')
                else:
                    f = text("TD",((startpos / aa2pixelratio) + 0.5 * ((endpos-startpos) / aa2pixelratio)), ((z * 60 ) + 4),fill='#2A0A29')
                if ((endpos-startpos) / aa2pixelratio) < 100 and ((endpos-startpos) / aa2pixelratio) >= 20:
                    myStyle.setFontSize('8px')
                    if "Thioesterase" in j:
                        f = text("TE",((startpos / aa2pixelratio) + 0.5 * ((endpos-startpos) / aa2pixelratio)), ((z * 60 ) + 3),fill='#2A0A29')
                    else:
                        f = text("TD",((startpos / aa2pixelratio) + 0.5 * ((endpos-startpos) / aa2pixelratio)), ((z * 60 ) + 4),fill='#2A0A29')
                elif ((endpos-startpos) / aa2pixelratio) < 20:
                    f = "notext"
                if f != "notext":
                    f.set_style(myStyle.getStyle())
                myStyle.setFontSize('12px')
                group = g()
                group.addElement(c)
                group.addElement(d)
                group.addElement(e)
                if f != "notext":
                    group.addElement(f)
                group.set_id("b" + str(qclusternr) + "_00%s"%w)
                s.addElement(group)
                tenr += 1
            else:
                c = LinearGradient("#929292","#DBDBDB","other_domain" + str(qclusternr) + "_" + str(othernr))
                d = LinearGradient("#DBDBDB","#929292","other_line" + str(qclusternr) + "_" + str(othernr))
                e = oh.createRect(str(10 + startpos / aa2pixelratio),str((z * 60 ) - 8),str((endpos-startpos) / aa2pixelratio),15,8,strokewidth=1,stroke='url(#other_line' + str(qclusternr) + "_" + str(othernr) + ")",fill="url(#other_domain" + str(qclusternr) + "_" + str(othernr) + ")")
                domname = (((((((((j.replace("0","")).replace("1","")).replace("2","")).replace("3","")).replace("4","")).replace("5","")).replace("6","")).replace("7","")).replace("8","")).replace("9","")
                if len(domname) == 1:
                    f = text(domname,((startpos / aa2pixelratio) + 0.5 * ((endpos-startpos) / aa2pixelratio)), ((z * 60 ) + 4),fill='#0B0B0B')
                elif len(domname) == 2:
                    f = text(domname,((-4 + startpos / aa2pixelratio) + 0.5 * ((endpos-startpos) / aa2pixelratio)), ((z * 60 ) + 4),fill='#0B0B0B')
                elif len(domname) == 3:
                    f = text(domname,((-12 + startpos / aa2pixelratio) + 0.5 * ((endpos-startpos) / aa2pixelratio)), ((z * 60 ) + 4),fill='#0B0B0B')
                if len(domname) > 3 or ((endpos-startpos) / aa2pixelratio) < 100:
                    myStyle.setFontSize('8px')
                    f = text(domname,((-16 + startpos / aa2pixelratio) + 0.5 * ((endpos-startpos) / aa2pixelratio)), ((z * 60 ) + 3),fill='#0B0B0B')
                if len(domname) > 4 and ((endpos-startpos) / aa2pixelratio) < 100:
                    myStyle.setFontSize('6px')
                    f = text(domname,((-16 + startpos / aa2pixelratio) + 0.5 * ((endpos-startpos) / aa2pixelratio)), ((z * 60 ) + 3),fill='#0B0B0B')
                if ((endpos-startpos) / aa2pixelratio) < 60:
                    f = "notext"
                if f != "notext":
                    f.set_style(myStyle.getStyle())
                myStyle.setFontSize('12px')
                group = g()
                group.addElement(c)
                group.addElement(d)
                group.addElement(e)
                if f != "notext":
                    group.addElement(f)
                group.set_id("b" + str(qclusternr) + "_00%s"%w)
                s.addElement(group)
                othernr += 1
            w += 1
        z += 1
    s.addElement(group)
    return s

def relativepositions(starts,ends,largestclustersize, screenwidth):
    rel_starts = []
    rel_ends = []
    #Assign relative start and end sites for visualization
    leftboundary = min([int(p) for p in starts+ends])
    for i in starts:
        i = float(float(int(i) - int(leftboundary)) / largestclustersize) * float(screenwidth * 0.75)
        i = int(i)
        rel_starts.append(i)
    for i in ends:
        i = float(float(int(i) - int(leftboundary)) / largestclustersize) * float(screenwidth * 0.75)
        i = int(i)
        rel_ends.append(i)
    return [rel_starts,rel_ends]


def geneclustersvg(genes,rel_starts,rel_ends,strands,geneposdict,pksnrpsprots,pksnrpsdomains,qclusternr, colors, screenwidth):
    nrgenes = len(genes)
    #Define relative start and end positions for plotting
    s = svg(x = 0, y = 0, width = (screenwidth * 0.75), height = (90 * len(pksnrpsprots)))
    viewbox = "0 -30 " + str(screenwidth * 0.8) + " " + str(60 * len(pksnrpsprots))
    s.set_viewBox(viewbox)
    s.set_preserveAspectRatio("none")

    #Add line behind gene arrows
    oh = ShapeBuilder()
    group = g()
    group.addElement(oh.createLine(10,60,10 + (screenwidth * 0.75),60, strokewidth = 2, stroke = "grey"))
    s.addElement(group)
    #Add gene arrows
    a = 0
    y = 0
    for x in range(nrgenes):
        group = g()
        #group.addElement(_gene_label(rel_starts[a],rel_ends[a],genes[a],y,screenwidth))
        group.addElement(_gene_arrow(10 + rel_starts[a],10 + rel_ends[a],strands[a],colors[a],65,10))
        #Can be used for domains
        #group.addElement(oh.createRect(rel_starts[a],45,(rel_ends[a]-rel_starts[a]),10, strokewidth = 2, stroke = "black", fill="#237845"))
        group.set_id("a" + str(qclusternr) + "_00%s"%x)
        s.addElement(group)
        if y == 0:
            y = 1
        elif y == 1:
            y = 0
        a += 1
    #Add NRPS/PKS domain depictions
    if len(pksnrpsprots) > 0:
        s = domains_svg_part(s, geneposdict,pksnrpsprots,pksnrpsdomains,qclusternr, screenwidth)
    return s

def create_clusterblast_svg(options, seq_record, screenwidth, searchtype="general"):
    #Create ClusterBlast svg
    clusterblastpositiondata = {}
    #Create alignment svg for each pair of hit&query
    geneclusters = utils.get_sorted_cluster_features(seq_record)
    if searchtype == "general":
        queryclusterdata = seq_record.queryclusterdata
        file_id = "clusterblast"
        logging.debug("Generating SVG for ClusterBlast")
    elif searchtype == "subclusters":
        queryclusterdata = seq_record.sc_queryclusterdata
        file_id = "subclusterblast"
        logging.debug("Generating SVG for SubClusterBlast")
    elif searchtype == "knownclusters":
        queryclusterdata = seq_record.kc_queryclusterdata
        file_id = "knownclusterblast"
        logging.debug("Generating SVG for KnownClusterBlast")
    idx = 1
    for genecluster in geneclusters:
        qclusternr = utils.get_cluster_number(genecluster)
#        print "Dumping queryclusterdata"
#        pprint.pprint(queryclusterdata)
        hitclusters = list(range(queryclusterdata[qclusternr][0] + 1))[1:]
        #Create svgs for pairwise gene cluster alignment
        colorschemedict,rgbcolorscheme = calculate_colorgroups(qclusternr, hitclusters, queryclusterdata, seq_record.internalhomologygroupsdict)
        for k in hitclusters:
            cresults = clusterblastresults(options, qclusternr, [k], queryclusterdata, colorschemedict, rgbcolorscheme, screenwidth, file_id, idx)
            s = cresults[0]
            clusterblastpositiondata[str(qclusternr) + "_"+str(k)] = cresults[1]
            outfile = open(options.svgdir + os.sep + file_id + str(qclusternr) + "_" + str(k) + ".svg","w")
            outfile.write(s.getXML())
            outfile.close()
        #Create svgs for multiple gene cluster alignment
        cresults = clusterblastresults(options, qclusternr, hitclusters, queryclusterdata, colorschemedict, rgbcolorscheme, screenwidth, file_id, idx)
        s = cresults[0]
        if s is not None:
            clusterblastpositiondata[str(qclusternr) + "_all"] = cresults[1]
            outfile = open(options.svgdir + os.sep + file_id + str(qclusternr) + "_all.svg","w")
            outfile.write(s.getXML())
            outfile.close()
        idx += 1
    return clusterblastpositiondata

def clusterblastresults(options, queryclusternumber,hitclusternumbers,queryclusterdata,colorschemedict,rgbcolorscheme, screenwidth, file_id, idx):
    #print "Generating svg for cluster",queryclusternumber
    #Extract data and generate color scheme
    nrhitclusters = queryclusterdata[queryclusternumber][0]
    hitclusterdata = queryclusterdata[queryclusternumber][1]
    if nrhitclusters == 0:
        return [None,[{},{},{}]]
    queryclustergenes = hitclusterdata[1][3]
    queryclustergenesdetails = hitclusterdata[1][4]
    colorgroupsdict = {}
    colorgroupslengthlist = []
    colorgroupslist = []
    for hitclusternumber in hitclusternumbers:
        colorgroups = hitclusterdata[hitclusternumber][0][hitclusternumber]
        colorgroupsdict[hitclusternumber] = colorgroups
        colorgroupslengthlist.append(len(colorgroups))
        colorgroupslist.append(colorgroups)
    #Find out whether hit gene cluster needs to be inverted compared to query gene cluster
    strandsbalancedict = findstrandsbalance(hitclusternumbers, hitclusterdata, queryclustergenes, queryclustergenesdetails, colorgroupsdict)
    #Calculate gene coordinates
    qstarts, qends, qstrands, qcolors, qnrgenes, hdata = generate_genecoords(strandsbalancedict, queryclustergenes, queryclustergenesdetails, colorschemedict, hitclusternumbers, hitclusterdata)
    #Find sizes of largest and smallest gene clusters of query & all hit clusters assessed
    largestclustersize, smallestclustersize, qclustersize = findclustersizes(options, hitclusternumbers, hdata, qstarts, qends)
    #Find relative positions
    qrelpositions = relativepositions(qstarts,qends,largestclustersize, screenwidth)
    qrel_starts = qrelpositions[0]
    qrel_ends = qrelpositions[1]
    qdata = [qrel_starts,qrel_ends,qstrands,qcolors]
    #Align gene clusters properly
    hdata, qdata = align_geneclusters(options, hitclusternumbers, hdata, qdata, largestclustersize, qclustersize, screenwidth)
    #Make SVG template
    width = 800
    height = 180 + (len(hitclusternumbers) * 50)
    s = svg(x=0, y=0, width=width, height=height)
    viewbox = "0 0 %d %d" % (width, height)
    s.set_viewBox(viewbox)
    s.set_preserveAspectRatio("none")
    #Construct gene arrows and lines which compose the gene cluster alignments
    s = clusterblast_gene_arrows(s, qnrgenes, qdata, hdata, rgbcolorscheme, queryclusternumber, hitclusternumbers, hitclusterdata, screenwidth, file_id, idx, strandsbalancedict)
    return [s,[qdata,hdata,strandsbalancedict]]

def calculate_colorgroups(queryclusternumber,hitclusternumbers,queryclusterdata,internalhomologygroupsdict):
    #Extract data and generate color scheme
    nrhitclusters = queryclusterdata[queryclusternumber][0]
    hitclusterdata = queryclusterdata[queryclusternumber][1]
    if nrhitclusters == 0:
        rgbcolorscheme = generate_rgbscheme(1)
        rgbcolorscheme.append("#FFFFFF")
        return {},rgbcolorscheme
    queryclustergenes = hitclusterdata[1][3]
    colorgroupsdict = {}
    colorgroupslengthlist = []
    colorgroupslist = []
    for hitclusternumber in hitclusternumbers:
        colorgroups = hitclusterdata[hitclusternumber][0][hitclusternumber]
        colorgroupsdict[hitclusternumber] = colorgroups
        colorgroupslengthlist.append(len(colorgroups))
        colorgroupslist.append(colorgroups)
    metacolorgroups = []
    internalgroups = internalhomologygroupsdict[queryclusternumber]
    for i in internalgroups:
        metagroup = []
        for j in i:
            for m in colorgroupslist:
                for l in m:
                    if j in l:
                        for k in l:
                            if k not in metagroup:
                                metagroup.append(k)
        if len(metagroup) > 1 and metagroup not in metacolorgroups:
            metacolorgroups.append(metagroup)
    #Generate RGB scheme
    rgbcolorscheme = generate_rgbscheme(len(metacolorgroups))
    rgbcolorscheme.append("#FFFFFF")
    #Create colorschemedict in which all genes that are hits of the same query gene get the same color
    colorschemedict = {}
    z = 0
    for i in queryclustergenes:
        for j in metacolorgroups:
            if i in j:
                for l in j:
                    if l in colorschemedict:
                        pass
                    else:
                        colorschemedict[l] = z
        if z in list(colorschemedict.values()):
            z += 1
    return colorschemedict,rgbcolorscheme

def findstrandsbalance(hitclusternumbers, hitclusterdata, queryclustergenes, queryclustergenesdetails, colorgroupsdict):
    #Finds out the average direction of an aligned gene cluster, to decide whether its orientation needs to be inverted.
    strandsbalancedict = {}
    for m in hitclusternumbers:
        hitclustergenesdetails = hitclusterdata[m][2]
        strandsbalance = 0
        for i in queryclustergenes:
            refstrand = queryclustergenesdetails[i][2]
            for j in colorgroupsdict[m]:
                if i in j:
                    for k in j:
                        if k in hitclusterdata[m][1] and hitclustergenesdetails[k][2] == refstrand:
                            strandsbalance += 1
                        elif k in hitclusterdata[m][1] and hitclusterdata[m][2][k][2] != refstrand:
                            strandsbalance = strandsbalance - 1
        strandsbalancedict[m] = strandsbalance
    return strandsbalancedict

def generate_genecoords(strandsbalancedict, queryclustergenes, queryclustergenesdetails, colorschemedict, hitclusternumbers, hitclusterdata):
    #Generate coordinates for SVG figure
    qnrgenes = len(queryclustergenes)
    qstarts =[]
    qends = []
    qstrands =[]
    qcolors = []
    for i in queryclustergenes:
        qgenedata = queryclustergenesdetails[i]
        if qgenedata[0] > qgenedata[1]:
            qstarts.append(qgenedata[0])
            qends.append(qgenedata[1])
        else:
            qstarts.append(qgenedata[1])
            qends.append(qgenedata[0])
        qstrands.append(qgenedata[2])
        if i in colorschemedict:
            qcolors.append(colorschemedict[i])
        else:
            qcolors.append("white")
    qstarts_ends = startendsitescheck(qstarts,qends)
    qstarts = qstarts_ends[0]
    qends = qstarts_ends[1]
    hdata = {}
    for m in hitclusternumbers:
        hitclustergenes = hitclusterdata[m][1]
        hitclustergenesdetails = hitclusterdata[m][2]
        hstarts =[]
        hends = []
        hstrands =[]
        hcolors = []
        for i in hitclustergenes:
            hgenedata = hitclustergenesdetails[i]
            if hgenedata[0] > hgenedata[1]:
                hstarts.append(hgenedata[0])
                hends.append(hgenedata[1])
            else:
                hstarts.append(hgenedata[1])
                hends.append(hgenedata[0])
            hstrands.append(hgenedata[2])
            if i in colorschemedict:
                hcolors.append(colorschemedict[i])
            else:
                hcolors.append("white")
        #Invert gene cluster if needed
        if strandsbalancedict[m] < 0:
            hstarts2 = []
            hends2 = []
            hstrands2 = []
            for i in hstarts:
                hstarts2.append(str(100000000 - int(i)))
            hstarts = hstarts2
            hstarts.reverse()
            for i in hends:
                hends2.append(str(100000000 - int(i)))
            hends = hends2
            hends.reverse()
            for i in hstrands:
                if i == "+":
                    hstrands2.append("-")
                elif i == "-":
                    hstrands2.append("+")
            hstrands = hstrands2
            hstrands.reverse()
            hcolors.reverse()
        hstarts_ends = startendsitescheck(hstarts,hends)
        hstarts = hstarts_ends[0]
        hends = hstarts_ends[1]
        hdata[m] = [hstarts,hends,hstrands,hcolors]
    return qstarts, qends, qstrands, qcolors, qnrgenes, hdata

def findclustersizes(options, hitclusternumbers, hdata, qstarts, qends):
    #Find cluster size of largest cluster
    clustersizes = []
    for m in hitclusternumbers:
        hclustersize = abs(int(hdata[m][1][-1]) - int(hdata[m][0][0]))
        clustersizes.append(hclustersize)
    qclustersize = abs(int(qends[-1]) - int(qstarts[0]))
    clustersizes.append(qclustersize)
    largestclustersize = max(clustersizes)
    if options.homologyscalelimit > 1.0:
        largestclustersize = min(largestclustersize, options.homologyscalelimit * qclustersize)
    smallestclustersize = min(clustersizes)
    return largestclustersize, smallestclustersize, qclustersize

def align_geneclusters(options, hitclusternumbers, hdata, qdata, largestclustersize, qclustersize, screenwidth):
    hdata2 = {}
    qdata2 = []
    qrel_starts = qdata[0]
    qrel_ends = qdata[1]
    query_modified = False
    for m in hitclusternumbers:
        hclustersize = abs(int(hdata[m][1][-1]) - int(hdata[m][0][0]))
        hrelpositions = relativepositions(hdata[m][0],hdata[m][1],largestclustersize, screenwidth)
        hrel_starts = hrelpositions[0]
        hrel_ends = hrelpositions[1]
        #Center-align smallest gene cluster
        if ( (options.homologyscalelimit > 1.0 and largestclustersize <= hclustersize ) or \
             (options.homologyscalelimit <= 1.0 and largestclustersize == hclustersize)) \
            and not query_modified:
            query_modified = True
            qrel_ends2 = []
            qrel_starts2 = []
            for i in qrel_starts:
                qrel_starts2.append(int(i) + int(float(float(float(largestclustersize - qclustersize) / 2.0) / largestclustersize) * float(screenwidth * 0.75)))
            for i in qrel_ends:
                qrel_ends2.append(int(i) + int(float(float(float(largestclustersize - qclustersize) / 2.0) / largestclustersize) * float(screenwidth * 0.75)))
            qrel_ends = qrel_ends2
            qrel_starts = qrel_starts2
        else:
            hrel_ends2 = []
            hrel_starts2 = []
            for i in hrel_starts:
                hrel_starts2.append(int(i) + int(float(float(float(largestclustersize - hclustersize) / 2.0) / largestclustersize) * float(screenwidth * 0.75)))
            for i in hrel_ends:
                hrel_ends2.append(int(i) + int(float(float(float(largestclustersize - hclustersize) / 2.0) / largestclustersize) * float(screenwidth * 0.75)))
            hrel_ends = hrel_ends2
            hrel_starts = hrel_starts2
        hdata2[m] = [hrel_starts,hrel_ends,hdata[m][2],hdata[m][3]]
        qdata2 = [qrel_starts,qrel_ends,qdata[2],qdata[3]]
    return hdata2, qdata2

def clusterblast_gene_arrows(s, qnrgenes, qdata, hdata, rgbcolorscheme, queryclusternumber, hitclusternumbers, hitclusterdata, screenwidth, file_id, idx, strandsbalancedict):
    queryclustergenes = hitclusterdata[1][3]
    queryclustergenesdetails = hitclusterdata[1][4]
    qrel_starts, qrel_ends, qstrands, qcolors = qdata
    #Add line behind query gene cluster gene arrows
    oh = ShapeBuilder()
    group = g()
    querycluster_acc = "Query sequence"
    acc = text(querycluster_acc, 5, 20)
    acc.set_class("clusterblast-acc")
    group.addElement(acc)
    group.addElement(oh.createLine(10,35,10 + (screenwidth * 0.75),35, strokewidth = 1, stroke = "grey"))
    group.setAttribute('label', querycluster_acc)
    group.set_class('clusterblast-cluster')
    s.addElement(group)
    #Add query gene cluster gene arrows
    a = 0
    for x in range(qnrgenes):
        #group.addElement(_gene_label(rel_starts[a],rel_ends[a],genes[a],y,screenwidth))
        if qcolors[a] == "white":
            arrow = _gene_arrow(10 + qrel_starts[a],10 + qrel_ends[a],qstrands[a],rgbcolorscheme[-1],40,10)
        else:
            arrow = _gene_arrow(10 + qrel_starts[a],10 + qrel_ends[a],qstrands[a],rgbcolorscheme[qcolors[a]],40,10)


        if len(hitclusternumbers) == 1:
            arrow.set_id("%s-%s_q%s_%s_%s" % (file_id, idx, queryclusternumber, hitclusternumbers[0], x))
        else:
            arrow.set_id("%s-%s_all_%s_0_%s" % (file_id, idx, queryclusternumber, x))
        arrow.set_class('clusterblast-orf')
        arrow.setAttribute('locus_tag', queryclustergenes[x])
        start, end, strand, desc = queryclustergenesdetails[queryclustergenes[x]]
        arrow.setAttribute('description', '%s[br]Location: %s - %s' % (desc, start, end))
        group.addElement(arrow)
        a += 1
    for m in hitclusternumbers:
        #Add line behind hit gene cluster gene arrows
        hitcluster_acc = hitclusterdata[m][6]
        hitcluster_desc = hitclusterdata[m][5][m-1].split('\t')[-1].replace('_', ' ')
        group = g()
        if hitcluster_acc.startswith('BGC'):
            acc = text('<a xlink:href="http://mibig.secondarymetabolites.org/repository/' + hitcluster_acc.partition("_c")[0] + '/index.html" target="_blank">'
                       + hitcluster_acc + '</a>: ' + hitcluster_desc, 5, 20 + 50 * (hitclusternumbers.index(m) + 1))
        else:
            acc = text(hitcluster_acc + ": " + hitcluster_desc, 5, 20 + 50 * (hitclusternumbers.index(m) + 1))
        acc.set_class("clusterblast-acc")
        group.addElement(acc)
       # desc = text(hitcluster_desc, (screenwidth * 0.45), 10 + 50 * (hitclusternumbers.index(m) + 1))
       # desc.set_class("clusterblast-desc")
       # group.addElement(desc)
        group.addElement(oh.createLine(10,35 + 50 * (hitclusternumbers.index(m) + 1),10 + (screenwidth * 0.75),35 + 50 * (hitclusternumbers.index(m) + 1), strokewidth = 1, stroke = "grey"))
        group.setAttribute('label', hitcluster_acc)
        group.setAttribute('description', hitcluster_desc)
        group.set_class('clusterblast-cluster')
        s.addElement(group)
        #Add hit gene cluster gene arrows
        hitclustergenes = hitclusterdata[m][1]
        hitclustergenesdetails = hitclusterdata[m][2]
        blastdetailsdict = hitclusterdata[m][-1]
        hits_accessions_dict = hitclusterdata[m][-2]
        if strandsbalancedict[m] < 0:
            hitclustergenes.reverse()
        hnrgenes = len(hitclustergenes)
        hrel_starts = hdata[m][0]
        hrel_ends = hdata[m][1]
        hstrands = hdata[m][2]
        hcolors = hdata[m][3]
        a = 0
        for x in range(hnrgenes):
            #group.addElement(_gene_label(rel_starts[a],rel_ends[a],genes[a],y,screenwidth))
            if hcolors[a] == "white":
                arrow = _gene_arrow(10 + hrel_starts[a],10 + hrel_ends[a],hstrands[a],rgbcolorscheme[-1],40 + 50 * (hitclusternumbers.index(m) + 1),10)
            else:
                arrow = _gene_arrow(10 + hrel_starts[a],10 + hrel_ends[a],hstrands[a],rgbcolorscheme[hcolors[a]],40 + 50 * (hitclusternumbers.index(m) + 1),10)
            if len(hitclusternumbers) == 1:
                arrow.set_id("%s-%s_h%s_%s_%s" % (file_id, idx, queryclusternumber, m, x))
            else:
                arrow.set_id("%s-%s_all_%s_%s_%s" % (file_id, idx, queryclusternumber, m, x))

            arrow.setAttribute('locus_tag', hitclustergenes[x])
            start, end, strand, desc, hitgene_acc = hitclustergenesdetails[hitclustergenes[x]]
            description_string = '%s[br]Location: %s - %s' % (desc, start, end)
            if hitclustergenes[x] in hits_accessions_dict:
                querygenes = hits_accessions_dict[hitclustergenes[x]]
                for querygene in querygenes:
                    evalue, bitscore, pid, pcov = blastdetailsdict[querygene + "_|_|_" + hitclustergenes[x]]
                    if pcov.isdigit() and int(pcov) > 100:
                        pcov = 100
                    pid = str(int(float(pid))) + "%"
                    pcov = str(int(float(pcov))) + "%"
                    description_string += '<br><br><b>BlastP hit with %s</b><br>Percentage identity: %s<br>Percentage coverage: %s<br>BLAST bit score: %s<br>E-value: %s' % (querygene, pid, pcov, bitscore, evalue)
            arrow.setAttribute('description', description_string)
            arrow.set_class("clusterblast-orf")
            group.addElement(arrow)
            a += 1
    return s
