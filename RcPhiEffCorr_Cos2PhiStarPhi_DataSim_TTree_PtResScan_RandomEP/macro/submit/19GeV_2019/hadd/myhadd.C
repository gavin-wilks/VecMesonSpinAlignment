#include <fstream>
#include <iostream>
#include <iomanip>

//----------------------------------------------------
void myhadd(TString list, TString jobId, TString output)
{
  TString command = 0;
  command.Append("hadd "+output+"_"+jobId+".root");
  //----------------------------------------------------------------------------------------------------
  // input
  if (list != NULL)   // if input file is ok
  {
    std::ifstream in;
    in.open(list,ios::in);  // input stream
    if(in)
    {
      cout << "input file list is ok" << endl;
      char str[255];       // char array for each file name
      while(in)
      {
	in.getline(str,255);  // take the lines of the file list
	if(str[0] != 0)
	{
	  TString addfile;
	  addfile = str;
          command.Append(" "+addfile);
	}
      }
    }
    std::cout << command.Data() << endl;
    gSystem->Exec(command.Data());
  }
}
