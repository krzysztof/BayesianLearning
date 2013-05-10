net Unnamed
{
 HEADER = 
  {
   ID = Unnamed;
   NAME = "Unnamed";
  };
 CREATION = 
  {
  };
 NUMSAMPLES = 1000;
 SCREEN = 
  {
   POSITION = 
    {
     CENTER_X = 0;
     CENTER_Y = 0;
     WIDTH = 76;
     HEIGHT = 36;
    };
   COLOR = 16250597;
   SELCOLOR = 12303291;
   FONT = 1;
   FONTCOLOR = 0;
   BORDERTHICKNESS = 3;
   BORDERCOLOR = 12255232;
  };
 WINDOWPOSITION = 
  {
   CENTER_X = 0;
   CENTER_Y = 0;
   WIDTH = 0;
   HEIGHT = 0;
  };
 BKCOLOR = 16777215;
 USER_PROPERTIES = 
  {
  };
 DOCUMENTATION = 
  {
  };
 SHOWAS = 3;

 node Smoker
  {
   TYPE = CPT;
   HEADER = 
    {
     ID = Smoker;
     NAME = "Smoker";
    };
   SCREEN = 
    {
     POSITION = 
      {
       CENTER_X = 0;
       CENTER_Y = 0;
       WIDTH = 76;
       HEIGHT = 36;
      };
     COLOR = 16250597;
     SELCOLOR = 12303291;
     FONT = 1;
     FONTCOLOR = 0;
     BORDERTHICKNESS = 1;
     BORDERCOLOR = 12255232;
    };
   USER_PROPERTIES = 
    {
    };
   DOCUMENTATION = 
    {
    };
   PARENTS = ();
   DEFINITION = 
    {
     NAMESTATES = (True, False);
     PROBABILITIES = (0.78800000, 0.21200000);
    };
   EXTRA_DEFINITION = 
    {
     DIAGNOSIS_TYPE = AUXILIARY;
     RANKED = FALSE;
     MANDATORY = FALSE;
     SETASDEFAULT = FALSE;
     SHOWAS = 4;
     FAULT_STATES = (0, 0);
     FAULT_NAMES = ("", "");
     FAULT_LABELS = ("", "");
     DEFAULT_STATE = 0;
     DOCUMENTATION = 
      {
      };
     DOCUMENTATION = 
      {
      };
     STATECOMMENTS = ("", "");
     STATEREPAIRINFO = ("", "");
     QUESTION = "";
    };
  };

 node Genetic
  {
   TYPE = CPT;
   HEADER = 
    {
     ID = Genetic;
     NAME = "Genetic";
    };
   SCREEN = 
    {
     POSITION = 
      {
       CENTER_X = 0;
       CENTER_Y = 0;
       WIDTH = 76;
       HEIGHT = 36;
      };
     COLOR = 16250597;
     SELCOLOR = 12303291;
     FONT = 1;
     FONTCOLOR = 0;
     BORDERTHICKNESS = 1;
     BORDERCOLOR = 12255232;
    };
   USER_PROPERTIES = 
    {
    };
   DOCUMENTATION = 
    {
    };
   PARENTS = ();
   DEFINITION = 
    {
     NAMESTATES = (True, False);
     PROBABILITIES = (0.70450000, 0.29550000);
    };
   EXTRA_DEFINITION = 
    {
     DIAGNOSIS_TYPE = AUXILIARY;
     RANKED = FALSE;
     MANDATORY = FALSE;
     SETASDEFAULT = FALSE;
     SHOWAS = 4;
     FAULT_STATES = (0, 0);
     FAULT_NAMES = ("", "");
     FAULT_LABELS = ("", "");
     DEFAULT_STATE = 0;
     DOCUMENTATION = 
      {
      };
     DOCUMENTATION = 
      {
      };
     STATECOMMENTS = ("", "");
     STATEREPAIRINFO = ("", "");
     QUESTION = "";
    };
  };

 node CoalWorker
  {
   TYPE = CPT;
   HEADER = 
    {
     ID = CoalWorker;
     NAME = "CoalWorker";
    };
   SCREEN = 
    {
     POSITION = 
      {
       CENTER_X = 0;
       CENTER_Y = 0;
       WIDTH = 76;
       HEIGHT = 36;
      };
     COLOR = 16250597;
     SELCOLOR = 12303291;
     FONT = 1;
     FONTCOLOR = 0;
     BORDERTHICKNESS = 1;
     BORDERCOLOR = 12255232;
    };
   USER_PROPERTIES = 
    {
    };
   DOCUMENTATION = 
    {
    };
   PARENTS = ();
   DEFINITION = 
    {
     NAMESTATES = (True, False);
     PROBABILITIES = (0.84550000, 0.15450000);
    };
   EXTRA_DEFINITION = 
    {
     DIAGNOSIS_TYPE = AUXILIARY;
     RANKED = FALSE;
     MANDATORY = FALSE;
     SETASDEFAULT = FALSE;
     SHOWAS = 4;
     FAULT_STATES = (0, 0);
     FAULT_NAMES = ("", "");
     FAULT_LABELS = ("", "");
     DEFAULT_STATE = 0;
     DOCUMENTATION = 
      {
      };
     DOCUMENTATION = 
      {
      };
     STATECOMMENTS = ("", "");
     STATEREPAIRINFO = ("", "");
     QUESTION = "";
    };
  };

 node BadDiet
  {
   TYPE = CPT;
   HEADER = 
    {
     ID = BadDiet;
     NAME = "BadDiet";
    };
   SCREEN = 
    {
     POSITION = 
      {
       CENTER_X = 0;
       CENTER_Y = 0;
       WIDTH = 76;
       HEIGHT = 36;
      };
     COLOR = 16250597;
     SELCOLOR = 12303291;
     FONT = 1;
     FONTCOLOR = 0;
     BORDERTHICKNESS = 1;
     BORDERCOLOR = 12255232;
    };
   USER_PROPERTIES = 
    {
    };
   DOCUMENTATION = 
    {
    };
   PARENTS = ();
   DEFINITION = 
    {
     NAMESTATES = (True, False);
     PROBABILITIES = (0.55050000, 0.44950000);
    };
   EXTRA_DEFINITION = 
    {
     DIAGNOSIS_TYPE = AUXILIARY;
     RANKED = FALSE;
     MANDATORY = FALSE;
     SETASDEFAULT = FALSE;
     SHOWAS = 4;
     FAULT_STATES = (0, 0);
     FAULT_NAMES = ("", "");
     FAULT_LABELS = ("", "");
     DEFAULT_STATE = 0;
     DOCUMENTATION = 
      {
      };
     DOCUMENTATION = 
      {
      };
     STATECOMMENTS = ("", "");
     STATEREPAIRINFO = ("", "");
     QUESTION = "";
    };
  };

 node LungCancer
  {
   TYPE = NOISY_MAX;
   HEADER = 
    {
     ID = LungCancer;
     NAME = "LungCancer";
    };
   SCREEN = 
    {
     POSITION = 
      {
       CENTER_X = 0;
       CENTER_Y = 0;
       WIDTH = 76;
       HEIGHT = 36;
      };
     COLOR = 16250597;
     SELCOLOR = 12303291;
     FONT = 1;
     FONTCOLOR = 0;
     BORDERTHICKNESS = 1;
     BORDERCOLOR = 12255232;
    };
   USER_PROPERTIES = 
    {
    };
   DOCUMENTATION = 
    {
    };
   PARENTS = (Smoker, Genetic, CoalWorker, BadDiet);
   DEFINITION = 
    {
     NAMESTATES = (True, False);
     STRENGTHS = (0, 1, 0, 1, 0, 1, 0, 1);
     PROBABILITIES = (0.53291268, 0.46708732, 0.00000000, 1.00000000, 
     0.29255624, 0.70744376, 0.00000000, 1.00000000, 0.17878671, 
     0.82121329, 0.00000000, 1.00000000, 0.06196542, 0.93803458, 
     0.00000000, 1.00000000, 0.00000000, 1.00000000);
    };
   EXTRA_DEFINITION = 
    {
     DIAGNOSIS_TYPE = AUXILIARY;
     RANKED = FALSE;
     MANDATORY = FALSE;
     SETASDEFAULT = FALSE;
     SHOWAS = 4;
     FAULT_STATES = (0, 0);
     FAULT_NAMES = ("", "");
     FAULT_LABELS = ("", "");
     DEFAULT_STATE = 0;
     DOCUMENTATION = 
      {
      };
     DOCUMENTATION = 
      {
      };
     STATECOMMENTS = ("", "");
     STATEREPAIRINFO = ("", "");
     QUESTION = "";
    };
  };
 OBSERVATION_COST = 
  {

   node Smoker
    {
     PARENTS = ();
     COSTS = (0.00000000);
    };

   node Genetic
    {
     PARENTS = ();
     COSTS = (0.00000000);
    };

   node CoalWorker
    {
     PARENTS = ();
     COSTS = (0.00000000);
    };

   node BadDiet
    {
     PARENTS = ();
     COSTS = (0.00000000);
    };

   node LungCancer
    {
     PARENTS = ();
     COSTS = (0.00000000);
    };
  };
};
