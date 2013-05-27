net Network1
{
 HEADER = 
  {
   ID = Network1;
   NAME = "Network1";
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

 node P1
  {
   TYPE = CPT;
   HEADER = 
    {
     ID = P1;
     NAME = "P1";
     COMMENT = "";
    };
   SCREEN = 
    {
     POSITION = 
      {
       CENTER_X = 107;
       CENTER_Y = 255;
       WIDTH = 83;
       HEIGHT = 58;
      };
     COLOR = 16250597;
     SELCOLOR = 12303291;
     FONT = 1;
     FONTCOLOR = 0;
     BORDERTHICKNESS = 1;
     BORDERCOLOR = 8388608;
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
     PROBABILITIES = (0.88481027, 0.11518973);
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

 node P2
  {
   TYPE = CPT;
   HEADER = 
    {
     ID = P2;
     NAME = "P2";
     COMMENT = "";
    };
   SCREEN = 
    {
     POSITION = 
      {
       CENTER_X = 317;
       CENTER_Y = 188;
       WIDTH = 91;
       HEIGHT = 70;
      };
     COLOR = 16250597;
     SELCOLOR = 12303291;
     FONT = 1;
     FONTCOLOR = 0;
     BORDERTHICKNESS = 1;
     BORDERCOLOR = 8388608;
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
     PROBABILITIES = (0.26228090, 0.73771910);
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

 node P3
  {
   TYPE = CPT;
   HEADER = 
    {
     ID = P3;
     NAME = "P3";
     COMMENT = "";
    };
   SCREEN = 
    {
     POSITION = 
      {
       CENTER_X = 488;
       CENTER_Y = 191;
       WIDTH = 69;
       HEIGHT = 43;
      };
     COLOR = 16250597;
     SELCOLOR = 12303291;
     FONT = 1;
     FONTCOLOR = 0;
     BORDERTHICKNESS = 1;
     BORDERCOLOR = 8388608;
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
     PROBABILITIES = (0.53210407, 0.46789593);
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

 node P4
  {
   TYPE = CPT;
   HEADER = 
    {
     ID = P4;
     NAME = "P4";
     COMMENT = "";
    };
   SCREEN = 
    {
     POSITION = 
      {
       CENTER_X = 686;
       CENTER_Y = 244;
       WIDTH = 52;
       HEIGHT = 32;
      };
     COLOR = 16250597;
     SELCOLOR = 12303291;
     FONT = 1;
     FONTCOLOR = 0;
     BORDERTHICKNESS = 1;
     BORDERCOLOR = 8388608;
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
     PROBABILITIES = (0.46669260, 0.53330740);
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

 node C1
  {
   TYPE = NOISY_MAX;
   HEADER = 
    {
     ID = C1;
     NAME = "C1";
     COMMENT = "";
    };
   SCREEN = 
    {
     POSITION = 
      {
       CENTER_X = 337;
       CENTER_Y = 399;
       WIDTH = 68;
       HEIGHT = 42;
      };
     COLOR = 16250597;
     SELCOLOR = 12303291;
     FONT = 1;
     FONTCOLOR = 0;
     BORDERTHICKNESS = 1;
     BORDERCOLOR = 8388608;
    };
   USER_PROPERTIES = 
    {
    };
   DOCUMENTATION = 
    {
    };
   PARENTS = (P1, P2, P3, P4);
   DEFINITION = 
    {
     NAMESTATES = (True, False);
     STRENGTHS = (0, 1, 0, 1, 0, 1, 0, 1);
     PROBABILITIES = (0.65573770, 0.34426230, 0.00000000, 1.00000000, 
     0.91666667, 0.08333333, 0.00000000, 1.00000000, 0.55555556, 
     0.44444444, 0.00000000, 1.00000000, 0.44444444, 0.55555556, 
     0.00000000, 1.00000000, 0.21276596, 0.78723404);
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

   node P1
    {
     PARENTS = ();
     COSTS = (0.00000000);
    };

   node P2
    {
     PARENTS = ();
     COSTS = (0.00000000);
    };

   node P3
    {
     PARENTS = ();
     COSTS = (0.00000000);
    };

   node P4
    {
     PARENTS = ();
     COSTS = (0.00000000);
    };

   node C1
    {
     PARENTS = ();
     COSTS = (0.00000000);
    };
  };
};
