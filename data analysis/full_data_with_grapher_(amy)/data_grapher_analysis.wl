(* ::Package:: *)

(* ::Input:: *)
(*(*graphing*)*)


(* ::Input:: *)
(*graphData[filename_,label_,color_]:=Module[{data=#[[1;;2]]&/@(Import[NotebookDirectory[]<>ToString[filename]])},ListPlot[data,Joined->True,PlotRange->{{0,2500},{-75,30}},PlotLabel->label,Frame->{{True,False},{False,False}},FrameLabel->{{"",None},{None,"vmax"}},FrameTicks->All,PlotStyle->color]]*)


(* ::Input:: *)
(*graphData["paired_data/paired-n0_1-k5","0.1 5", RGBColor[0,0,0]] *)


(* ::Input:: *)
(*#[[1;;2]]&/@(Import[NotebookDirectory[]<>ToString["paired_data/paired-n0_1-k5.txt"]])*)


(* ::Input:: *)
(*(*graphing two data sets on one set of axes*)*)


(* ::Input:: *)
(*data1=#[[1;;2]]&/@(Import[NotebookDirectory[]<>ToString["pnkcc_data/pnkcc_0"]]);*)
(*data2=#[[1;;2]]&/@(Import[NotebookDirectory[]<>ToString["pnkcc_data/pnkcc_05"]]);*)


(* ::Input:: *)
(*ListPlot[{data1,data2},Joined->True,PlotRange->{{0,2500},{-75,30}}]*)


(* ::Input:: *)
(*(*graph kcc2*)*)
(*inhibited= RGBColor[61/255,143/255,186/255]*)
(*default=RGBColor[100/255,100/255,100/255]*)
(*excited = RGBColor[255/255,81/255,133/255]*)
(*graph=GraphicsGrid[{*)
(*{graphData["vmax_data/vmax_0_01","0.01\[Cross]KCC2",excited]},{graphData["vmax_data/vmax_0_5","0.5\[Cross]KCC2", excited]},*)
(*{graphData["vmax_data/vmax_1","1\[Cross]KCC2", default]},{graphData["vmax_data/vmax_1_5","1.5\[Cross]KCC2", inhibited]}, *)
(*{graphData["vmax_data/vmax_2","2\[Cross]KCC2", inhibited]}, {graphData["vmax_data/vmax_5","5\[Cross]KCC2",inhibited]}*)
(*}]*)
(*Export["transparent.png",graph,Background->None]*)


(* ::Input:: *)
(*SystemOpen["transparent.png"]*)


(* ::Input:: *)
(*(*graph nkcc1*)*)
(*nkcc1graph = GraphicsGrid[{*)
(*{graphData["pnkcc_data/pnkcc_2","2\[Cross]NKCC1", excited]},*)
(*{graphData["pnkcc_data/pnkcc_1_5","1.5\[Cross]NKCC1", excited]},*)
(*{graphData["pnkcc_data/pnkcc_1","1\[Cross]NKCC1", default]},*)
(*{graphData["pnkcc_data/pnkcc_075","0.75\[Cross]NKCC1", inhibited]}, *)
(*{graphData["pnkcc_data/pnkcc_05","0.5\[Cross]NKCC1", inhibited]}, *)
(*{graphData["pnkcc_data/pnkcc_01","0.1\[Cross]NKCC1", inhibited]}*)
(*}]*)
(*Export["transparent.png",nkcc1graph,Background->None]*)


(* ::Input:: *)
(*SystemOpen["transparent.png"]*)


(* ::Input:: *)
(*SystemOpen["transparent.png"]*)


(* ::Input:: *)
(*SystemOpen["transparent.png"]*)


(* ::Input:: *)
(**)


(* ::Input:: *)
(**)


(* ::Input:: *)
(*(*graph paired*)*)
(*GraphicsGrid[{*)
(*{graphData["paired_data/paired-n0_1-k5","0.1 5", Blue]},*)
(*{graphData["paired_data/paired-n0_25-k2","0.25 2", Blue]},*)
(*(*0.4 1.6 -> no spikes*)*)
(*{graphData["paired_data/paired-n0_5-k1_5","0.5 1.5", Blue]},*)
(*{graphData["paired_data/paired-n0_6-k1_4","0.6 1.4", Blue]},*)
(*(*0.8 1.2 -> no spikes*)*)
(*(*0.9 1.1 -> ? *)*)
(*{graphData["paired_data/paired-n1-k1","1 1", Blue]},*)
(*(*1.1 0.9*)*)
(*{graphData["paired_data/paired-n1_2-k0_8", "1.2 0.8", Red]},*)
(*{graphData["paired_data/paired-n1_5-k0_5","1.5 0.5", Blue]},*)
(*{graphData["paired_data/paired-n2-k0_25","2 0.25", Blue]},*)
(*{graphData["paired_data/paired-n5-k0_1","5 0.1", Blue]}*)
(*}]*)


(* ::Input:: *)
(**)
