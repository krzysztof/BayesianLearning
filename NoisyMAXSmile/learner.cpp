#include <iostream>
#include <fstream>
#include "smile/smile.h"
#include "smile/smilearn.h"
#include <stdlib.h>
#include <string>
#include <vector>
#include <time.h>

using namespace std;

void printCPT(DSL_node *node)
{
   DSL_network* net = node->Network(); // node network                                                                   
   int handle = node->Handle();
   DSL_nodeDefinition *def = node->Definition();
   const DSL_Dmatrix &cpt = *def->GetMatrix();
   const DSL_idArray &outcomes = *def->GetOutcomesNames();
   const DSL_intArray &parents = net->GetParents(handle);
   int parentCount = parents.NumItems();

   DSL_intArray coords;
   for (int elemIdx = 0; elemIdx < cpt.GetSize(); elemIdx ++)
   {
	  string name = "";
      cpt.IndexToCoordinates(elemIdx, coords);
      //cout << "P(" << node->GetId() << " = " << outcomes[coords[parentCount]] << " | ";
      for (int parentIdx = 0; parentIdx < parentCount; parentIdx ++)
      {
         DSL_node *parentNode = net->GetNode(parents[parentIdx]);
         const DSL_idArray &parentStates = *parentNode->Definition()->GetOutcomesNames();
         //cout << parentNode->GetId() << " = " << parentStates[coords[parentIdx]];             

		 name += parentStates[coords[parentIdx]];
		 name += " ";
      }
	  name += outcomes[coords[parentCount]];
	  cout << name << " " << cpt[elemIdx] << endl;
      //cout << ") = " << cpt[elemIdx] << endl;
   }
}

DSL_network LearnParams(string data_infile, string network_infile){
  DSL_dataset ds;
  if (ds.ReadFile(data_infile.c_str()) != DSL_OKAY) {
     cout << "Cannot read data file... exiting." << endl;
     exit(1);
  }

    DSL_network net;
  if (net.ReadFile(network_infile.c_str(), DSL_XDSL_FORMAT) != DSL_OKAY) {
     cout << "Cannot read network... exiting." << endl;
     exit(1);
  }

  vector<DSL_datasetMatch> matches;
  string err;
  if (ds.MatchNetwork(net, matches, err) != DSL_OKAY) {
     cout << "Cannot match network... exiting." << endl;
     exit(1);
  }

  DSL_em em;
  em.SetUniformizeParameters(true);
  if (em.Learn(ds, net, matches) != DSL_OKAY) {
     cout << "Cannot learn parameters... exiting." << endl;
     exit(1);
  }
  //net.WriteFile("res_tut_5.xdsl", DSL_XDSL_FORMAT);
  return net;

}

int main(int argc, char* argv[])
{
	string data_infile = argv[1];
	string network_infile = argv[2];
	DSL_network net = LearnParams(data_infile, network_infile);
	int cancer = net.FindNode("LungCancer");
	DSL_node* node = net.GetNode(cancer);
	if(argv[3] == string("EM")) {
		printCPT(node);
	} 
	else if(argv[3] == string("Genie")) {
		net.GetNode(cancer)->ChangeType(DSL_NOISY_MAX);
		printCPT(node);
	} else if(argv[3] == string("Naive")) {
		// naiwne uczenie
	} else if(argv[3] == string("Gauss")) {
		// gauss-jordan
	}

	return 0;
}

