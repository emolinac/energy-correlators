# energy-correlations
## Code Details
- In the DecayTree ntuples, the jets are organizes in events. Each event has a number of candidates and, almost always, the first candidate in the event
is the hardest jet. Therefore, at the beginning of each jet loop, there is a block that makes sure to only check the first candidate in the event:
```C++
mcrecotree->GetEntry(evt);
if (evt != 0)
{
    if (mcrecotree->eventNumber != last_eventNum) maxjetpT_found = false;
    if (last_eventNum == mcrecotree->eventNumber) continue;
}
last_eventNum = mcrecotree->eventNumber;
if (maxjetpT_found) continue;
```
