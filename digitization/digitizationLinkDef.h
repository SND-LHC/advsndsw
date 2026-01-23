#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ nestedclasses;
#pragma link C++ nestedtypedef;
#pragma link C++ namespace advsnd;
#pragma link C++ defined_in namespace advsnd;

#pragma link C++ class advsnd::DigitizePoints<AdvTargetPoint, AdvTargetHit>+;
#pragma link C++ class advsnd::DigitizePoints<AdvMuFilterPoint, AdvMuFilterHit>+;
#pragma link C++ class advsnd::LinkPointsToDigi<AdvTargetPoint>+;
#pragma link C++ class advsnd::LinkPointsToDigi<AdvMuFilterPoint>+;

#pragma link C++ function advsnd::GetStripId<AdvTargetPoint>;
#pragma link C++ function advsnd::GetStripId<AdvMuFilterPoint>;
#pragma link C++ function advsnd::Digitize<AdvTargetPoint, AdvTargetHit>;
#pragma link C++ function advsnd::Digitize<AdvMuFilterPoint, AdvMuFilterHit>;
#pragma link C++ function advsnd::McLink<AdvTargetPoint>;
#pragma link C++ function advsnd::McLink<AdvMuFilterPoint>;

#endif