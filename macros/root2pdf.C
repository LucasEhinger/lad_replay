#include <TFile.h>
#include <TCanvas.h>
#include <TDirectory.h>
#include <TH1.h>
#include <TKey.h>
#include <TString.h>
#include <iostream>
#include <functional>

void root2pdf(const char* rootFileName, const char* pdfFileName = "", int n = 2, int m = 2) {
  // If pdfFileName is empty, derive it from rootFileName
  TString pdfFile(pdfFileName);
  if (pdfFile.IsNull()) {
    pdfFile = TString(rootFileName);
    pdfFile.ReplaceAll(".root", ".pdf");
    pdfFileName = pdfFile.Data();
  }

  // Open the ROOT file
  TFile* file = TFile::Open(rootFileName);
  if (!file || file->IsZombie()) {
    std::cerr << "Error: Cannot open file " << rootFileName << std::endl;
    return;
  }

  // Create a canvas for drawing
  TCanvas canvas("canvas", "Canvas", 800, 600);
  canvas.Divide(m, n); // Divide canvas into m columns and n rows

  // Start the PDF
  canvas.Print((TString(pdfFileName) + "[").Data());

  // Recursive function to loop through directories
  std::function<void(TDirectory*)> processDirectory = [&](TDirectory* dir) {
    TList* keys = dir->GetListOfKeys();
    if (!keys) return;

    int pad = 1;

    for (int i = 0; i < keys->GetSize(); ++i) {
      TKey* key = (TKey*)keys->At(i);
      if (!key) continue;

      TObject* obj = key->ReadObj();
      if (obj->InheritsFrom(TH1::Class())) {
        // If it's a histogram, draw it
        canvas.cd(pad);
        obj->Draw();
        pad++;

        // If the canvas is full, print the page and reset
        if (pad > n * m) {
          canvas.Print(pdfFileName);
          canvas.Clear();
          canvas.Divide(m, n);
          pad = 1;
        }
      } else if (obj->InheritsFrom(TCanvas::Class())) {
        // If it's a TCanvas, print it directly
        TCanvas* objCanvas = (TCanvas*)obj;
        objCanvas->Print(pdfFileName);
      } else if (obj->InheritsFrom(TDirectory::Class())) {
        // If it's a directory, process it recursively
        processDirectory((TDirectory*)obj);
      }
    }

    // Print remaining histograms on the last page
    if (pad > 1) {
      canvas.Print(pdfFileName);
      canvas.Clear();
      canvas.Divide(m, n);
    }
  };

  // Start processing from the root directory
  processDirectory(file);

  // Close the PDF
  canvas.Print((TString(pdfFileName) + "]").Data());

  // Close the ROOT file
  file->Close();
  delete file;

  std::cout << "Histograms and canvases saved to " << pdfFileName << std::endl;
}