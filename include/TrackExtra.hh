#ifndef TrackExtra_h
#define TrackExtra_h 1

#include "globals.hh"
#include "G4ThreeVector.hh"
#include "G4ParticleDefinition.hh"
#include "G4Track.hh"
#include "G4Allocator.hh"
#include "G4VUserTrackInformation.hh"

class TrackExtra : public G4VUserTrackInformation
{
public:
    TrackExtra();
    TrackExtra(const G4Track* aTrack);
    TrackExtra(const TrackExtra* aTrackInfo);
    virtual ~TrackExtra();
    
    inline void *operator new(size_t);
    inline void operator delete(void *aTrackInfo);
    inline int operator ==(const TrackExtra& right) const
    {return (this==&right);}
    
    void Print() const;
    
private:
    G4int                 originalTrackID;
    
public:
    inline G4int GetOriginalTrackID() const {return originalTrackID;}
};

extern G4Allocator<TrackExtra> aTrackExtraAllocator;

inline void* TrackExtra::operator new(size_t)
{ void* aTrackInfo;
    aTrackInfo = (void*)aTrackExtraAllocator.MallocSingle();
    return aTrackInfo;
}

inline void TrackExtra::operator delete(void *aTrackInfo)
{ aTrackExtraAllocator.FreeSingle((TrackExtra*)aTrackInfo);}

#endif
