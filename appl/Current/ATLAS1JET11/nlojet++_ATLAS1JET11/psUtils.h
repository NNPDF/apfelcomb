
TGraph* psMakeBand(THStack* inputStack, TString option = "quadrature")
{
  //
  //  makes graph with errors out of stack of histograms
  //

  int debug = 1;
  if (!inputStack)
    {
      cout <<"histogram stack is not found "<< endl;
      return 0;
    }
  
  TList*  localList = inputStack->GetHists();
  const Int_t listSize = localList->GetSize();

  TGraphErrors* centralValue = (TGraphErrors*)psTH12TGraph((TH1D*)localList->At(0));
  TGraphAsymmErrors* graphBand = new TGraphAsymmErrors();


  Int_t switchCase = 0;
  Double_t clFactor = 1.64485;
  Double_t bandFactor = 1.;
  option.ToLower();

  if ( option.Contains("abkm") || option.Contains("alek") ) 
    {
      switchCase = 50;
      bandFactor = clFactor;
    }
  if (option.Contains("gjr")) 
    {
      switchCase = 10;
      bandFactor = clFactor;
    }
  if (option.Contains("nnpdf")) 
    {
      switchCase = 100;
      bandFactor = clFactor;
    }
  if (option.Contains("hera")) 
    {
      switchCase = 10;
      bandFactor = clFactor;
    }
  if (option.Contains("envelope")) switchCase = 1111;
  if (option.Contains("cteq") || option.Contains("ct10")) switchCase = 10;
  if (option.Contains("mrst") || option.Contains("mstw")) switchCase = 10;
  if (option.Contains("quadrature")) switchCase = 50;
  
  for (int iPoint = 1; iPoint < centralValue->GetN(); iPoint++)
    {
      Double_t otherYValues[1000];
      Double_t centralXValue, centralXValueError, centralYValue, centralYValueErrorL, centralYValueErrorH;
      centralValue->GetPoint(iPoint, centralXValue, centralYValue);
      centralXValueError = 0.5*centralValue->GetErrorX(iPoint);

      graphBand->SetPoint(iPoint, centralXValue, centralYValue);
      graphBand->SetPointEXhigh(iPoint, centralXValueError);
      graphBand->SetPointEXlow (iPoint, centralXValueError);

      for (int iHist = 1; iHist < listSize; iHist++)
	{
	  TH1D* htemp = (TH1D*)localList->At(iHist);
	  otherYValues[iHist] = htemp->GetBinContent(htemp->FindBin(centralXValue));
	}

      switch (switchCase)
	{
	case (1111):
	  calculateErrorsEnvelope(listSize, otherYValues, 
				  centralYValue, 
				  centralYValueErrorL, 
				  centralYValueErrorH);

	  graphBand->SetPointEYhigh(iPoint, centralYValueErrorH);
	  graphBand->SetPointEYlow (iPoint, centralYValueErrorL);
	  break;

	case (100):
	      calculateErrorsNNPDF(listSize, otherYValues, centralYValue, centralYValueErrorH);
	      graphBand->SetPointEYhigh(iPoint, clFactor * centralYValueErrorH);
	      graphBand->SetPointEYlow (iPoint, clFactor * centralYValueErrorH);
	  break;
	case (10):
	  if (option.Contains("asym"))
	    {
	      calculateErrorsHessianAsymmetric(listSize, otherYValues, 
					       centralYValue, 
					       centralYValueErrorL, 
					       centralYValueErrorH);
	      graphBand->SetPointEYhigh(iPoint, bandFactor * centralYValueErrorH);
	      graphBand->SetPointEYlow (iPoint, bandFactor * centralYValueErrorL);
	    }
	  else
	    {
	      calculateErrorsHessianSymmetric(listSize, otherYValues, centralYValueErrorH);
	      graphBand->SetPointEYhigh(iPoint, bandFactor * centralYValueErrorH);
	      graphBand->SetPointEYlow (iPoint, bandFactor * centralYValueErrorH);
	    }
	  break;
	  // 	case (20):
	  // 	  break;
	  
	  // 	case (30):
	  // 	  break;
	case (50):
	  calculateErrorsSymmetric(listSize, otherYValues, centralYValue, centralYValueErrorH);
	  graphBand->SetPointEYhigh(iPoint, centralYValueErrorH);
	  graphBand->SetPointEYlow (iPoint, centralYValueErrorH);

	  break;
	default:
	  //calculateErrorsSymmetric(listSize, otherYValues, centralYValue, centralYValueErrorH);
	  centralYValueErrorH = 2. *centralYValue;  // just to test
	  graphBand->SetPointEYhigh(iPoint, centralYValueErrorH);
	  graphBand->SetPointEYlow (iPoint, centralYValueErrorH);
	}

    }
  //graphBand->Print();
  //  cout<<"\tmax = "<<graphBand->GetMaximum()<<"\t";
  //  cout<<"\tmin = "<<graphBand->GetMinimum()<<endl<<endl;
  return graphBand;
}

void  calculateErrorsEnvelope(Int_t nValues, Double_t *oValues, Double_t cValue, Double_t &cValueErrorL, Double_t &cValueErrorH)
{
  
  cValueErrorL = 0.;
  cValueErrorH = 0.;
  
  for (int i = 1; i < nValues ; i++)
    {
      //cout <<"i = "<<i<<" cValue = "<<cValue<<" var = "<<oValues[i]<< endl;
      Double_t diff = oValues[i] - cValue;
      if (diff >= 0.)
	  cValueErrorH = TMath::Max(cValueErrorH, diff);
      else
	  cValueErrorL = TMath::Max(cValueErrorL, -diff);
    }
  //cout <<"errorL = "<<cValueErrorL<<" errorH = "<<cValueErrorH<< endl;
  return;
}
void  calculateErrorsHessianAsymmetric(Int_t nValues, Double_t *oValues, Double_t cValue, Double_t &cValueErrorL, Double_t &cValueErrorH)
{

  cValueErrorL = 0.;
  cValueErrorH = 0.;

  for (int i = 1; i < nValues ; i += 2)
    {
      Double_t diffPlus = oValues[i] - cValue;
      Double_t diffMinus = oValues[i+1] - cValue;

      Double_t currentErrorPlus = 0.;
      currentErrorPlus = TMath::Max(currentErrorPlus, diffPlus);
      currentErrorPlus = TMath::Max(currentErrorPlus, diffMinus);

      Double_t currentErrorMinus = 0.;
      currentErrorMinus = TMath::Max(currentErrorMinus, -diffPlus);
      currentErrorMinus = TMath::Max(currentErrorMinus, -diffMinus);

      cValueErrorL += currentErrorMinus * currentErrorMinus;
      cValueErrorH += currentErrorPlus * currentErrorPlus;
    }

  cValueErrorL = TMath::Sqrt(cValueErrorL);
  cValueErrorH = TMath::Sqrt(cValueErrorH);
}

void  calculateErrorsHessianSymmetric(Int_t nValues, Double_t *oValues, Double_t &cValueError)
{

  cValueError = 0.;

  for (int i = 1; i < nValues ; i += 2)
    {
      Double_t diff = oValues[i] - oValues[i+1];

      cValueError += diff * diff;
    }

  cValueError = 0.5 * TMath::Sqrt(cValueError);
}


void calculateErrorsSymmetric(Int_t nValues, Double_t *oValues, Double_t cValue, Double_t &cValueError)
{
  cValueError = 0.;

  for (int i = 1; i < nValues ; i++)
    {
      Double_t diff = cValue - oValues[i];
      cValueError += diff*diff;
    }

  cValueError = 0.5 * std::sqrt(cValueError);
}

void calculateErrorsAsymmetric(Int_t nValues, Double_t *oValues, Double_t cValue, Double_t &cValueErrorL, Double_t &cValueErrorH)
{
  cValueErrorL = 0.;
  cValueErrorH = 0.;

  for (int i = 1; i < nValues ; i++)
    {
      Double_t diff = cValue - oValues[i];
      if (diff >= 0.)
	{
	  cValueErrorH += diff*diff;
	}
      else
	{
	  cValueErrorL += diff*diff;
	}
    }

  cValueErrorL = TMath::Sqrt(cValueErrorL);
  cValueErrorH = TMath::Sqrt(cValueErrorH);
}

void calculateErrorsNNPDF(Int_t nValues, Double_t *oValues, Double_t cValue, Double_t &cValueError)
{
  cValueError = 0.;
  for (int i = 1; i < nValues ; i++)
    {
      Double_t diff = cValue - oValues[i];
      cValueError += diff*diff;
    }

  cValueError = TMath::Sqrt(cValueError/(nValues - 2.));

}


TGraph* psTH12TGraph(TH1 *h1, bool binError = false)
{
  //
  // convert the histogram h1 into a graph
  //
  if (!h1) 
    {
      cout <<name<< " histogram not found !" << endl;
      return 0;
    }
  
  TGraphErrors* g1 = new TGraphErrors();
  
  Double_t x, y, ex, ey;
  for (Int_t i = 1; i <= h1->GetNbinsX(); i++) 
    {
      x =  h1->GetBinCenter(i);
      ex = h1->GetBinWidth(i);
      //--
      y = h1->GetBinContent(i);
      (binError) ? ey = h1->GetBinError(i): ey = 0;
      //--
      g1->SetPoint(i,x,y);
      g1->SetPointError(i,ex,ey);
    }
  return g1;
}
//
TGraphAsymmErrors* psAddBands(TGraphAsymmErrors* g1, TGraphAsymmErrors* g2)
{
  TGraphAsymmErrors* totalBand = new TGraphAsymmErrors();
  if (g1->GetN() != g2->GetN())
    {
      cout<<"!!! ERROR in addBands(gr*, gr*)"<<endl;
      cout <<" number of points in graph #1 "<<g1->GetN()<<"is not equal to nember of points in graph #2 "<<g2->GetN()<<endl;
      return totalBand;
    }
  
  for (int iPoint = 1; iPoint < g1->GetN(); iPoint++)
    {
      Double_t xValue, yValue;
      Double_t xValueErrorL, xValueErrorH;
      Double_t yValueErrorL, yValueErrorH, yValueErrorL1, yValueErrorH1, yValueErrorL2, yValueErrorH2;
      
      g1->GetPoint(iPoint, xValue, yValue);
      xValueErrorL = g1->GetErrorXlow(iPoint);
      xValueErrorH = g1->GetErrorXhigh(iPoint);
      
      yValueErrorL1 = g1->GetErrorYlow(iPoint);
      yValueErrorH1 = g1->GetErrorYhigh(iPoint);
      yValueErrorL2 = g2->GetErrorYlow(iPoint);
      yValueErrorH2 = g2->GetErrorYhigh(iPoint);
      
      yValueErrorL = TMath::Sqrt(yValueErrorL1 * yValueErrorL1 + yValueErrorL2 * yValueErrorL2);
      yValueErrorH = TMath::Sqrt(yValueErrorH1 * yValueErrorH1 + yValueErrorH2 * yValueErrorH2);
      
      totalBand->SetPoint(iPoint, xValue, yValue);
      totalBand->SetPointEXhigh(iPoint, xValueErrorH);
      totalBand->SetPointEXlow (iPoint, xValueErrorL);
      totalBand->SetPointEYhigh(iPoint, yValueErrorH);
      totalBand->SetPointEYlow (iPoint, yValueErrorL);
    }

  return totalBand;
}
  //
void getBandRange(TGraphAsymmErrors* g, Double_t &MIN, Double_t &MAX)
{
  MIN = 1000.0;
  MAX = -1000.0;
  for (Int_t iPoint = 1; iPoint < g->GetN(); iPoint++)
    {
      Double_t xValue, yValue;
      Double_t yValueErrorHigh, yValueErrorLow;
      g->GetPoint(iPoint, xValue, yValue);
      yValueErrorHigh = g->GetErrorYhigh(iPoint);
      yValueErrorLow  = g->GetErrorYlow (iPoint);
      MIN = TMath::Min(MIN, yValue - yValueErrorLow);
      MAX = TMath::Max(MAX, yValue + yValueErrorHigh);
    }

}

//
//
//
