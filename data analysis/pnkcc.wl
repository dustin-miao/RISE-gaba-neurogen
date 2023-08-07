(* ::Package:: *)

(* ::Input:: *)
(*pnkccVals={0,0.01,0.025,0.1,0.2,0.3,0.4,0.5,0.6,0.75,0.8,1,1.2,1.3,1.5,1.6,1.8,2};*)


(* ::Input:: *)
(*data2d=#[[All,1;;2]]&/@(Import[NotebookDirectory[]<>"pnkcc/"<>ToString[#]]&/@pnkccVals);*)
(*data3d=MapIndexed[Insert[#1,pnkccVals[[#2[[1]]]],1]&,data2d,{2}];*)
(*points=Flatten[data3d,1];*)


(* ::Input:: *)
(*countSpikesFile[filename_,threshold_:0,dtime_:6]:=Module[*)
(*{data=#[[1;;2]]&/@(Import[NotebookDirectory[]<>ToString[filename]]),count=0,prev=-dtime},*)
(*Do[If[data[[i,2]]>threshold&&data[[i,1]]-prev>dtime,prev=data[[i,1]];count+=1],{i,1,Length[data]}];*)
(*count*)
(*]*)


(* ::Input:: *)
(*countSpikesList[data_,threshold_:0,dtime_:6]:=Module[*)
(*{count=0,prev=-dtime},*)
(*Do[If[data[[i,3]]>threshold&&data[[i,2]]-prev>dtime,prev=data[[i,2]];count+=1],{i,1,Length[data]}];*)
(*count*)
(*]*)


(* ::Input:: *)
(*numSpikes=countSpikesList[#]&/@data3d*)
(*maxSpikes=Max[numSpikes]*)


(* ::Input:: *)
(*ListLogLinearPlot[Transpose[{pnkccVals,countSpikesFile["pnkcc/"<>ToString[#]]&/@pnkccVals}],ScalingFunctions->{"None","Log"},Joined->True,AxesLabel->{"pnkcc","Spikes"},PlotMarkers->{Automatic}]*)


(* ::Input:: *)
(*pnkcc2=ListPlot[Transpose[{pnkccVals,countSpikesFile["pnkcc/"<>ToString[#]]&/@pnkccVals}],Joined->True,PlotMarkers->{Automatic},PlotStyle->{Thickness[0.0125],inhibitedColor}]*)


(* ::Input:: *)
(*inhibitedColor=RGBColor[61/255,143/255,186/255];*)
(*defaultColor=RGBColor[100/255,100/255,100/255];*)
(*excitedColor=RGBColor[255/255,81/255,133/255];*)
(*colorFunctionNew[x_]:=Blend[{inhibitedColor,excitedColor},x]*)


(* ::Input:: *)
(*graph=Show[#]&@*)
(*MapIndexed[ListPointPlot3D[#1,*)
(*ColorFunction->Function[{x,y,z},colorFunctionNew[numSpikes[[#2[[1]]]]/maxSpikes]],*)
(*PlotRange->{{0,2},{0,2500},{-90,30}},*)
(*PlotStyle->{PointSize[0.005]},*)
(*ScalingFunctions->{None,"Reverse",None},*)
(*AxesLabel->{"pnkcc","time","vm"},*)
(*ViewPoint->{-Pi,-Pi,-Pi}*)
(*]&,*)
(*data3d,*)
(*{1}*)
(*]*)


(* ::Input:: *)


(* ::Input:: *)
(*Export["pnkcc.png",pnkcc2,Background->None]*)


(* ::Input:: *)
(*graph*)