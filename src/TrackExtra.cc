#include "TrackExtra.hh"
#include "G4ios.hh"

G4Allocator<TrackExtra> aTrackExtraAllocator;

TrackExtra::TrackExtra()
{
    originalTrackID = 0;
}

TrackExtra::TrackExtra(const G4Track* aTrack)
{
    originalTrackID = aTrack->GetTrackID();
}

TrackExtra::TrackExtra(const TrackExtra* aTrackInfo)
{
    originalTrackID = aTrackInfo->originalTrackID;
}

TrackExtra::~TrackExtra(){;}

void TrackExtra::Print() const
{
    
}
