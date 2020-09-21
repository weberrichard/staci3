#include "Graph.h"

//--------------------------------------------------------------
vector<vector<int> > segmenter(vector<int> ev)
{
  vector<vector<int> > everySegment;
  vector<int> segment;

  while(ev.size() !=0)
  {
    segment.push_back(ev[0]);
    segment.push_back(ev[1]);
    ev.erase(ev.begin());
    ev.erase(ev.begin());

    for(int j=0; j<segment.size(); j+=2)
    {
      for(int i=0; i<ev.size(); i+=2)
      {
        if(ev[i] == segment[j] || ev[i+1] == segment[j] || ev[i] == segment[j+1] || ev[i+1] == segment[j+1])
        {
          segment.push_back(ev[i]);
          segment.push_back(ev[i+1]);
          ev.erase(ev.begin() + i);
          ev.erase(ev.begin() + i);
          j=-2; // setting back j to the first item that is 0th
          break;
        }
      }
    }
    everySegment.push_back(segment);
    segment.clear();
  }

  return everySegment;
}
