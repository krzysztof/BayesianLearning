net Network4
{
 HEADER = 
  {
   ID = Network4;
   NAME = "Network4";
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
     PROBABILITIES = (0.33000000, 0.67000000);
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
       CENTER_X = 246;
       CENTER_Y = 216;
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
     PROBABILITIES = (0.22000000, 0.78000000);
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
       CENTER_X = 390;
       CENTER_Y = 212;
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
     PROBABILITIES = (0.01000000, 0.99000000);
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
       CENTER_X = 527;
       CENTER_Y = 222;
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
     PROBABILITIES = (0.55000000, 0.45000000);
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

 node P5
  {
   TYPE = CPT;
   HEADER = 
    {
     ID = P5;
     NAME = "P5";
     COMMENT = "";
    };
   SCREEN = 
    {
     POSITION = 
      {
       CENTER_X = 665;
       CENTER_Y = 267;
       WIDTH = 33;
       HEIGHT = 20;
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
     PROBABILITIES = (0.30000000, 0.70000000);
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

 node P6
  {
   TYPE = CPT;
   HEADER = 
    {
     ID = P6;
     NAME = "P6";
     COMMENT = "";
    };
   SCREEN = 
    {
     POSITION = 
      {
       CENTER_X = 761;
       CENTER_Y = 363;
       WIDTH = 33;
       HEIGHT = 20;
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
     PROBABILITIES = (0.66000000, 0.34000000);
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
       CENTER_X = 388;
       CENTER_Y = 431;
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
   PARENTS = (P1, P2, P3, P4, P5, P6);
   DEFINITION = 
    {
     NAMESTATES = (True, False);
     STRENGTHS = (0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1);
     PROBABILITIES = (0.46355653, 0.53644347, 0.00000000, 1.00000000, 
     0.12486160, 0.87513840, 0.00000000, 1.00000000, 0.00000000, 
     6.18367347, 0.00000000, 1.00000000, 0.11430932, 0.88569068, 
     0.00000000, 1.00000000, 0.16502206, 0.83497794, 0.00000000, 
     1.00000000, 0.90982726, 0.09017274, 0.00000000, 1.00000000, 
     0.00000000, 1.07846012);
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

   node P5
    {
     PARENTS = ();
     COSTS = (0.00000000);
    };

   node P6
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
