#include <iostream>
#include <fstream>
#include "smile/smile.h"
#include <stdlib.h>
#include <string>
#include <vector>
#include <time.h>
#include <map>

using namespace std;

typedef map<string, string> keyMap;
typedef map<string, double> valMap;
typedef map< keyMap, valMap > cptMap;

void printMap(keyMap &keymap) {
    keyMap::iterator it;
    for (it = keymap.begin(); it!= keymap.end() ; ++it) {
        cout << ":"<<(*it).first << ": = :" << (*it).second <<":" << endl;
    }
}

cptMap get_cptmap(DSL_node *node) {
    DSL_network* net = node->Network(); // node network                                                                   
    int handle = node->Handle();
    DSL_nodeDefinition *def = node->Definition();
    const DSL_Dmatrix &cpt = *def->GetMatrix();
    const DSL_idArray &outcomes = *def->GetOutcomesNames();
    const DSL_intArray &parents = net->GetParents(handle);
    int parentCount = parents.NumItems();

    DSL_intArray coords;
    cptMap cptmap;
    for (int elemIdx = 0; elemIdx < cpt.GetSize(); elemIdx ++) {
        cpt.IndexToCoordinates(elemIdx, coords);
        //cout << "P(" << node->GetId() << " = " << outcomes[coords[parentCount]] << " | ";
        keyMap km;
        for (int parentIdx = 0; parentIdx < parentCount; parentIdx ++) {
         //   if (parentIdx > 0) cout << ", ";
            DSL_node *parentNode = net->GetNode(parents[parentIdx]);
            const DSL_idArray &parentStates = *parentNode->Definition()->GetOutcomesNames();
            km[string(parentNode->GetId())] = string(parentStates[coords[parentIdx]]);
            //cout << parentNode->GetId() << " = " << parentStates[coords[parentIdx]];             
        }
        //printMap(km);
        if (cptmap.count(km) == 0){
            valMap vm;
            vm[string(outcomes[coords[parentCount]])] = cpt[elemIdx];
            cptmap[km] = vm;
        } else {
            cptmap[km][string(outcomes[coords[parentCount]])] = cpt[elemIdx];
        }

        //cout << ") = " << cpt[elemIdx] << endl;
    }

    return cptmap;
}

void generateDatafile(int n, string inname) {
	DSL_network theNet;
	theNet.ReadFile(inname.c_str());
	int cancer = theNet.FindNode("LungCancer");
	DSL_node *node = theNet.GetNode(cancer);
    cptMap cpt1 = get_cptmap(node);
    keyMap km1;
    km1[string("Smoker")] = string("True");
    km1[string("Genetic")] = string("True");
    km1[string("CoalWorker")] = string("True");
    km1[string("BadDiet")] = string("True");
    printMap(km1);
    cout << "SIZE: "<< cpt1.size() << endl;
    cout << cpt1[km1]["True"];
}

int main(int argc, char* argv[])
{
   //cancerNet();
   //generateDatafile(100, "CancerOR.xdsl", "CancerOR_generatedbyCpp.txt");
    cptMap supermap;
    valMap vm1;
    keyMap km1;
    km1[string("AA")] = string("aa");
    km1[string("BB")] = string("bb");

    keyMap km2;
    km2[string("BB")] = string("bb");
    km2[string("AA")] = string("aa");

    vm1[string("F")] = 0.667;
    vm1[string("T")] = 0.333;
    supermap[km1] = vm1;

    cout << supermap[km2]["T"] << endl;

    srand(time(NULL));
    generateDatafile(atoi(argv[1]), argv[2]);
    return 0;
}
