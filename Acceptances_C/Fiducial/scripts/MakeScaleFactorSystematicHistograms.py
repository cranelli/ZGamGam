#Python Code for looking at the Histograms in a root
#file, and printing the content of each bin.

import sys
import Config
from ROOT import TFile
from ROOT import TH1F
#from decimal import getcontext

#HIST_DIR="../Histograms/LepGammaGammaFinalElandMuUnblindAll_2015_5_27_ScaleFactor_PU_PDF_Reweights/"


CHANNELS=Config.CHANNELS
HIST_TYPES=Config.HIST_TYPES
ACCEPTANCE_TYPES=Config.ACCEPTANCE_TYPES
SFS=Config.SFS
DIRS = ["UP", "DN"]

def MakeScaleFactorSystematicHistograms():
    file_loc = Config.ACCEPTANCE_FILE_LOC
    file = TFile(file_loc, 'READ')

    outfile_loc = Config.SCALEFACTOR_SYSTEMATICS_FILE_LOC
    outfile=TFile(outfile_loc, 'RECREATE')
    
    #list = file.GetListOfKeys()
    for channel in CHANNELS:
        for acceptance_type in ACCEPTANCE_TYPES:
            for dir in DIRS:
                for hist_type in HIST_TYPES:
                    MakeIndividualSFSystematicHistograms(file, outfile, channel, acceptance_type, dir, hist_type)
                    MakeTotalSFSystematicHistograms(file, outfile, channel, acceptance_type, dir, hist_type)
            
# Goes over each channel's scale factors, and makes the ScaleFactor Syst Histogram.
# Just the difference between the two histograms
def MakeIndividualSFSystematicHistograms(file, outfile, channel, acceptance_type, dir, hist_type):
    expectedHistName=channel+"_Acceptance_"+acceptance_type+"_"+hist_type
    expectedHist = file.Get(expectedHistName)
    for scalefactor in SFS[channel]:
        # Histogram of Acceptances with raised, lopwered ScaleFactor
        scalefactorHistName = channel+"_"+scalefactor+dir+"_Acceptance_"+acceptance_type+"_"+hist_type
        print scalefactorHistName
        scalefactorHist = file.Get(scalefactorHistName)
            
        scalefactorSystematicHist = scalefactorHist.Clone(channel+"_"+acceptance_type+"_"+scalefactor+"_ScaleFactorSystematic"+dir+"_"+hist_type)
        scalefactorSystematicHist.Reset("ICE")
        # Loop over All bins in the Histogram, including the overflow bin.
        for bin_index in range(1, expectedHist.GetNbinsX()+2):
            scalefactorSyst = CalcAbsDifference(scalefactorHist, expectedHist, bin_index)
            scalefactorSystematicHist.SetBinContent(bin_index, scalefactorSyst)
                
        scalefactorSystematicHist.Print()
        scalefactorSystematicHist.Write()

# Adds the systematics for each scale factor in quadrature to get the total scale factor systematic.
# Stores in a Histogram
def MakeTotalSFSystematicHistograms(file, outfile, channel, acceptance_type, dir, hist_type):
    expectedHistName=channel+"_Acceptance_"+acceptance_type+"_"+hist_type
    expectedHist = file.Get(expectedHistName)

    scalefactorTotalSystematicHistName=channel+"_"+acceptance_type+"_ScaleFactorTotalSystematic"+dir+"_"+hist_type
    scalefactorTotalSystematicHist = expectedHist.Clone(scalefactorTotalSystematicHistName)
    scalefactorTotalSystematicHist.Reset("ICE")

    # Include Overflow
    for bin_index in range(1, scalefactorTotalSystematicHist.GetNbinsX()+2):
        differences=[]
        for scalefactor in SFS[channel]:
            # Histogram of Acceptances with raised, lowered ScaleFactor
            scalefactorHistName = channel+"_"+scalefactor+dir+"_Acceptance_"+acceptance_type+"_"+hist_type
            scalefactorHist = file.Get(scalefactorHistName)
            differences.append(CalcAbsDifference(scalefactorHist, expectedHist, bin_index))
            
        scalefactorTotalSyst = sumInQuadrature(differences)
        scalefactorTotalSystematicHist.SetBinContent(bin_index, scalefactorTotalSyst)

    scalefactorTotalSystematicHist.Print()
    scalefactorTotalSystematicHist.Write()


# Take the difference between the bins' of two histograms, and returns the absolute of the difference.
def CalcAbsDifference(hist1, hist2, bin_index):
    difference = hist1.GetBinContent(bin_index)-hist2.GetBinContent(bin_index)
    return abs(difference)

#Sum the differences in Quadrature
def sumInQuadrature(differences):
    sumdif2=0
    for difference in differences:
        sumdif2 += difference **2
    rootsumdif2= sumdif2 **0.5
    return rootsumdif2

"""
scalefactorSystUPHistName=channel+"_"+scalefactor+"UP"+acceptance_type+"_scalefactorSystematicUP_"+histType
        scalefactorSystUPHist=centralPDFHist.Clone(scalefactorSystUPHistName)
        scalefactorSystUPHist.Reset("ICE")
    
        scalefactorSystDNHistName=channel+"_"+acceptance_type+"_scalefactorSystematicDN_"+histType
        scalefactorSystDNHist=centralPDFHist.Clone(scalefactorSystDNHistName)
        scalefactorSystDNHist.Reset("ICE")




    # Select Histogram with Expected Values
    
    print expectedHist.GetName()
    outfileUP =open(outFileDir+"SystErrorsUP_"+expectedHist.GetName()+".dat", 'w')
    outfileDN =open(outFileDir+"SystErrorsDN_"+expectedHist.GetName()+".dat", 'w')
    # Store the Difference Systematics Together in a Dictionairy
            differencesUP=dict()
            differencesDN=dict()
            # Loop Over Histograms with modified Scale Factor Values
            for key in list:
                histName = key.GetName()
                if channel+"_"+histType in histName:
                    print histName
                    hist = file.Get(histName)
                    if hist.ClassName() == "TH1F":
                        # Histogram Names
                        printh1BinDifference(hist, expectedHist)
                    elif hist.ClassName() == "TH2F":
                        #printh2BinDifference(hist, expectedHist)
                        if "UP" in histName:
                            differencesUP[histName]=h2CalcDifferences(hist, expectedHist)
                        elif "DN" in histName:
                            differencesDN[histName]=h2CalcDifferences(hist, expectedHist)
                    else:
                        print "This type of Class Not Supported for Printing Bin Content"
            #Differences Up
            outputUP=""
            for key in differencesUP:
                print key
            # print differencesUP[key]
            systematicsUP= sumInQuadrature(differencesUP)
            for systematicUP in systematicsUP:
                outputUP+=str(format(systematicUP, '.6f')).rjust(10)
            outfileUP.write(outputUP)
            # print sumInQuadrature(differencesUP)
            
            # Differences Down
            outputDN=""
            for key in differencesDN:
                print key
            #    print differencesDN[key]
            systematicsDN= sumInQuadrature(differencesDN)
            for systematicDN in systematicsDN:
                outputDN+=str(format(systematicDN, '.6f')).rjust(10)
            outfileDN.write(outputDN)
            # print sumInQuadrature(differencesDN)
            outfileUP.close()
            outfileDN.close()
    file.Close()

# Sum the differences, for each list of differences in the Dictionairy. 
def sumInQuadrature(differencesDict):
    sumdif2=[]
    for key in differencesDict:
        differences = differencesDict[key]
        for i in range(len(differences)):
            if len(sumdif2) != len(differences) :
                sumdif2.insert(i, differences[i] **2)
            else:
                sumdif2[i] += differences[i] **2
    #Than take the square root of each element
    rootsumdif2=[]
    for i in range(len(sumdif2)):
        rootsumdif2.insert(i, sumdif2[i] **0.5)
    return rootsumdif2
    

# For a two 1D Histograms, prints the differences of the bin contents
def printh1BinDifference(h1, h1Expected):
    #hist = TH1F()
    output=""
    for i in range(1, h1.GetNbinsX()+1):
        output = str(h1.GetBinLowEdge(i)) +" : " +str(h1.GetBinContent(i)-h1Expected.GetBinContent(i))
        #print output
        #print h1.GetBinContent(i)

def h2CalcDifferences(h2, h2Expected):
    differences=[]
    for y in range(1, h2.GetNbinsY()+1):
        # Loop Over Columns (Photons' Location)
        for x in range(1, h2.GetNbinsX()+1):
            differences.append(h2.GetBinContent(x,y)-h2Expected.GetBinContent(x,y))
    return differences
        
"""
# def printh2BinDifference(h2, h2Expected):
#outfile =open(outFileDir+"BinDifferences_"+h2.GetName()+".html", 'w')
#print "Bins Content:"
#outfile.write(htmlTableHeader)
#output=""
#htmlOutput=""
    
# Loop Over Rows (Pt)
#for j in range(1, h2.GetNbinsY()+1):
#outfile.write("""<tr align="center">""")
# Loop Over Columns (Photons' Location)
#for i in range(1, h2.GetNbinsX()+1):
#output = str(h2.GetXaxis().GetBinLowEdge(i))+","+str(h2.GetYaxis().GetBinLowEdge(j))+" : "+ str(format(h2.GetBinContent(i,j),'.3f'))+" pm " + str(format(h2.GetBinError(i,j),'.3f'))
            #print output
            # Bin Content and Errors in HTML format
#htmlOutput = """<td valign="middle">""" + str(format(h2.GetBinContent(i,j)-h2Expected.GetBinContent(i,j),'.4f'))+ """</td>"""
#outfile.write(htmlOutput)
#outfile.write("""</tr>""")
            
#outfile.write(htmlTableCloser)
#outfile.close()
    
                                                                                                 
    
if __name__=="__main__":
    MakeScaleFactorSystematicHistograms()
