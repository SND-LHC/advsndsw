#include "AdvHit.h"

#include <iostream>

// -----   Public method Print   -------------------------------------------
void AdvHit::Print() const
{
    std::cout << " AdvHit: in detector " << detector_id_ << "\tsignal: " << signal_ << " \ttime: " << time_ << "\n";
}
