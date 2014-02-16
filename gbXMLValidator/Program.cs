using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Xml;
using System.IO;
using DOEgbXML;
using VectorMath;

namespace gbXMLValidator
{
    class Program
    {
        //define tolerances for the tests
        static double coordtol = DOEgbXMLBasics.Tolerances.coordToleranceIP;
        //path of Text file output
        static string outputpath = @"C:\gbXML";
        //a good time to figure out what to do with my units

        //an example main that has all of the tests that we would like to perform
        static void Main(string[] args)
        {

            //Basic File Preparation and Checks --------------------------------------------------------------------------------------------
            
            DOEgbXML.XMLParser parser = new XMLParser();
            //1-load the file
            //2-check for valid XML and valid gbXML against the XSD

            //hardcoded for testing purposes, eventually this file will be uploaded, or sent via a Restful API call
            string path = "C:\\Users\\Chiensi\\Documents\\C\\CarmelSoft\\gbXML Project\\Phase 2\\Test Files\\Test Case 5 - Standard File.xml";
            XmlReader xmlreader = XmlReader.Create(path);
            XmlDocument myxml = new XmlDocument();
            myxml.Load(xmlreader);

            //3-get the namespace
            XmlNamespaceManager nsm = parser.getnsmanager(myxml);
            //figure out if metric or USIP (we have not found a reason to use this yet)
            parser.getunits(nsm, myxml);




            //Begin Parsing the XML and reporting on it----------------------------------------------------------
            //make a reporting object
            DOEgbXMLReportingObj report = new DOEgbXMLReportingObj();

            //Basic Uniqueness Constraints check-------------------------------------------------

            //ensure that all names of spaces are unique
            report = DOEgbXML.gbXMLSpaces.UniqueSpaceIdTest(myxml, nsm, report);
            //process report
            report.Clear();
            //ensure that all space boundary names are unique
            report = DOEgbXML.gbXMLSpaces.UniqueSpaceBoundaryIdTest(myxml, nsm, report);
            //process report
            report.Clear();
            
            //Unique CAD Object IDs?
            
            //Space Tests
            //3-check for non-planar objects for all Spaces' polyloops
            
            
            List<DOEgbXML.gbXMLSpaces> spaces = DOEgbXML.gbXMLSpaces.getSimpleSpaces(myxml, nsm);

            report = DOEgbXML.gbXMLSpaces.SpaceSurfacesPlanarTest(spaces, report);
            //process report
            report.Clear();

            //4-check for self-intersecting polygons
            
            //Vertex Matching------------------------------------------------

            //try to parse out the Space Boundary polyloops
            Dictionary<Vector.CartCoord,Tuple<List<string>,List<bool>>> sbvertices = GetSGVertices(nsm, myxml);
            //do all vertices have at least one match?  If yes, PASS and move on, if not, then.

            VertexListToFile(sbvertices,"ShellGeometryCoords.txt");
            //try to parse out the surfaces into my surface objects
            //check the vertex files
            Dictionary<Vector.CartCoord, Tuple<List<string>, List<bool>>> surfvertices = GetSurfVertices(nsm, myxml);
            VertexListToFile(surfvertices, "SurfacesCoords.txt");
            //check the vertex files
            report = gbXMLSpaces.findStraySBVertices(@"C:\Temp\gbXML\SurfacesCoords.txt", report);
            //process report
            report.Clear();

            //Surface tests----------------------------------------------------------------------------------
            //Basic Requirements ------------------------------------------------------

            //Are there at least 4 surface definitions?  (see the surface requirements at the campus node)
            report = SurfaceDefinitions.AtLeast4Surfaces(myxml, nsm, report);
            //process report
            report.Clear();
            //Does the AdjacentSpaceId not exceed the max number allowable?
            report = SurfaceDefinitions.AtMost2SpaceAdjId(myxml, nsm, report);
            //process report
            report.Clear();
            //Are all required elements and attributes in place?
            report = SurfaceDefinitions.RequiredSurfaceFields(myxml, nsm, report);
            //process report
            report.Clear();
            //ensure that all names of surfaces are unique
            report = DOEgbXML.SurfaceDefinitions.SurfaceIDUniquenessTest(myxml, nsm, report);
            //process report
            report.Clear();
            //Does the id of the surface match the id of its corresponding space boundary?
            //Does the polyloop of the surface match the polyloop of the space boundary?
            //need to determine whether this is a metric or an IP coordinate!!
            report.tolerance = DOEgbXMLBasics.Tolerances.coordToleranceIP;
            report = SurfaceDefinitions.SurfaceMatchesSpaceBoundary(myxml, nsm, report);
            //process report
            report.Clear();

            //Does the polyloop right hand rule vector form the proper azimuth and tilt? (with and without a CADModelAzimuth)
            report.tolerance = DOEgbXMLBasics.Tolerances.VectorAngleTolerance;

            //Is the Lower Left Corner properly defined?

            //planar Surface Check
            List<SurfaceDefinitions> surfaces = DOEgbXML.XMLParser.MakeSurfaceList(myxml, nsm);
            SurfaceDefinitions.TestSurfacePlanarTest(surfaces);

            //self-intersecting Polygon Check
            report.testSummary = "This test checks to ensure the polygon definition for the surface forms non-intersecting enclosed polygon";
            List<Vector.MemorySafe_CartCoord> coordlist = new List<Vector.MemorySafe_CartCoord>();
            foreach (SurfaceDefinitions surface in surfaces)
            {
                for (int i = 0; i < surface.PlCoords.Count; i++)
                {
                    Vector.MemorySafe_CartCoord c = new Vector.MemorySafe_CartCoord(surface.PlCoords[i].X,surface.PlCoords[i].Y,surface.PlCoords[i].Z);
                    coordlist.Add(c);
                }
                bool r = Vector.BruteForceIntersectionTest(coordlist);
                if (!r)
                {
                    report.TestPassedDict.Add(surface.SurfaceId, false);
                }
            }

            //Vertex Matching


            //Openings Tests-----------------------------------------------------

            //Shading Devices Tests----------------------------------------------


            

        }

        private static Dictionary<Vector.CartCoord, Tuple<List<string>, List<bool>>> GetSGVertices(XmlNamespaceManager nsm, XmlDocument zexml)
        {
            Dictionary<Vector.CartCoord,Tuple<List<string>,List<bool>>> sbdict = new Dictionary<Vector.CartCoord, Tuple<List<string>,List<bool>>>();

            //get each space boundary
            XmlNodeList nodes = zexml.SelectNodes("/gbXMLv5:gbXML/gbXMLv5:Campus/gbXMLv5:Building/gbXMLv5:Space/gbXMLv5:SpaceBoundary", nsm);

            int sbcount = 0;
            foreach (XmlNode sbnode in nodes)
            {
                //record the name of the space boundary
                string sbname="";
                XmlAttributeCollection spaceAtts = sbnode.Attributes;
                foreach (XmlAttribute at in spaceAtts)
                {
                    if (at.Name == "surfaceIdRef")
                    {
                        sbname = at.Value;
                        break;
                    }
                }
                //create list of unique vertices
                XmlNode closedshell = sbnode.FirstChild;
                XmlNode PL = closedshell.FirstChild;

                foreach (XmlNode cp in PL)
                {
                    
                    if (cp.Name == "CartesianPoint")
                    {
                        sbdict = AddCoordinate(sbcount, sbname, sbdict, cp);
                    }
                }
                //XmlNodeList cnodes = sbnode.SelectNodes("/gbXMLv5:gbXML/gbXMLv5:Campus/gbXMLv5:Building/gbXMLv5:Space/gbXMLv5:SpaceBoundary/gbXMLv5:PlanarGeometry//gbXMLv5:PolyLoop", nsm);
                sbcount++;
            }
            return sbdict;
        }

        private static Dictionary<Vector.CartCoord, Tuple<List<string>, List<bool>>> GetSurfVertices(XmlNamespaceManager nsm, XmlDocument zexml)
        {
            Dictionary<Vector.CartCoord, Tuple<List<string>, List<bool>>> surfdict = new Dictionary<Vector.CartCoord, Tuple<List<string>, List<bool>>>();

            //get each space boundary
            XmlNodeList nodes = zexml.SelectNodes("/gbXMLv5:gbXML/gbXMLv5:Campus/gbXMLv5:Surface", nsm);

            int surfcount = 0;
            foreach (XmlNode surfnode in nodes)
            {
                //record the name of the space boundary
                string surfname = "";
                XmlAttributeCollection spaceAtts = surfnode.Attributes;
                foreach (XmlAttribute at in spaceAtts)
                {
                    if (at.Name == "id")
                    {
                        surfname = at.Value;
                        break;
                    }
                }
                //create list of unique vertices
                foreach (XmlNode childnode in surfnode)
                {
                    if (childnode.Name == "PlanarGeometry")
                    {
                        XmlNode PL = childnode.FirstChild;
                        foreach (XmlNode cp in PL)
                        {

                            if (cp.Name == "CartesianPoint")
                            {
                                surfdict = AddCoordinate(surfcount, surfname, surfdict, cp);
                            }
                        }
                        surfcount++;
                    }
                }
            }
            return surfdict;
        }

        private static Vector.CartCoord makeCoord(XmlNode cartesianPoint)
        {
            Vector.CartCoord coordinst = new Vector.CartCoord();
            if (cartesianPoint.HasChildNodes)
            {
                XmlNodeList coordinates = cartesianPoint.ChildNodes;
                int pointCount = 1;
                foreach (XmlNode coordinate in coordinates)
                {
                    switch (pointCount)
                    {
                        case 1:
                            coordinst.X = Convert.ToDouble(coordinate.InnerText);
                            break;
                        case 2:
                            coordinst.Y = Convert.ToDouble(coordinate.InnerText);
                            break;
                        case 3:
                            coordinst.Z = Convert.ToDouble(coordinate.InnerText);
                            break;
                    }
                    pointCount++;
                }
            }
            return coordinst;
        }

        private static bool isCoordUnique(Dictionary<Vector.CartCoord,Tuple<List<string>,List<bool>>> clist, Vector.CartCoord coord)
        {
            foreach (Vector.CartCoord approvedcoord in clist.Keys)
            {
                //a text for tolerances?
                if (approvedcoord.X == coord.X && approvedcoord.Y == coord.Y && approvedcoord.Z == coord.Z)
                {
                    return false;
                }
            }
            return true;
        }

        private static Dictionary<Vector.CartCoord, Tuple<List<string>,List<bool>>> 
            AddCoordinate(int surfacecount, string sbname, Dictionary<Vector.CartCoord, Tuple<List<string>, List<bool>>> surfdict, XmlNode cp)
        {
            Vector.CartCoord retcoord = makeCoord(cp);
            if (surfacecount == 0)
            {
                if (isCoordUnique(surfdict, retcoord))
                {
                    List<string> sbnames = new List<string>();
                    sbnames.Add(sbname);
                    List<bool> bl = new List<bool>();
                    bl.Add(true);
                    var t = Tuple.Create(sbnames, bl);
                    surfdict.Add(retcoord, t);
                }

            }
            else
            {
                if (isCoordUnique(surfdict, retcoord))
                {
                    List<string> sbnames = new List<string>();
                    sbnames.Add(sbname);
                    List<bool> bl = new List<bool>();
                    bl.Add(true);
                    var t = Tuple.Create(sbnames, bl);
                    surfdict.Add(retcoord, t);
                }
                else
                {
                    //this could be MAXIMIZED a little bit by activating the continue currently disabled below
                    foreach (KeyValuePair<Vector.CartCoord, Tuple<List<string>, List<bool>>> kp in surfdict)
                    {
                        //temporary
                        double dx = Math.Abs(retcoord.X - kp.Key.X);
                        double dy = Math.Abs(retcoord.Y - kp.Key.Y);
                        double dz = Math.Abs(retcoord.Z - kp.Key.Z);

                        if (dx <= coordtol && dy <= coordtol && dz <= coordtol)
                        {
                            if (dx == 0 && dy == 0 && dz == 0)
                            {
                                surfdict[kp.Key].Item1.Add(sbname);
                                surfdict[kp.Key].Item2.Add(true);

                            }
                            else
                            {
                                surfdict[kp.Key].Item1.Add(sbname);
                                surfdict[kp.Key].Item2.Add(false);
                            }
                            //continue;
                        }
                        else
                        {

                        }
                    }
                }
            }
            return surfdict;
        }

        private static void VertexListToFile(Dictionary<Vector.CartCoord, Tuple<List<string>, List<bool>>> vertices, string filename)
        {
            string filepath = Path.Combine(outputpath, Path.GetFileName(filename));
            using (TextWriter tw = File.CreateText(filepath))
            {
                string header = "Coordinate;Shared Surfaces Array;Perfect Match Array";
                tw.WriteLine(header);
                foreach (KeyValuePair<Vector.CartCoord, Tuple<List<string>, List<bool>>> kp in vertices)
                {
                    string s = "[";
                    s += kp.Key.X.ToString() + ",";
                    s += kp.Key.Y.ToString() + ",";
                    s += kp.Key.Z.ToString() + "];";

                    int surfcount = 0;
                    foreach (string surface in kp.Value.Item1)
                    {
                        if (surfcount == 0)
                        {
                            s += "[";
                            s += surface + ",";
                        }
                        else if (surfcount == kp.Value.Item1.Count - 1)
                        {
                            s += surface + "];";
                        }
                        else
                        {
                            s += surface + ",";
                        }
                        surfcount++;

                    }
                    int blcount = 0;
                    foreach (bool bl in kp.Value.Item2)
                    {
                        if (blcount == 0)
                        {
                            s += "[";
                            s += bl.ToString() + ",";
                        }
                        else if (blcount == kp.Value.Item2.Count - 1)
                        {
                            s += bl.ToString() + "];";
                        }
                        else
                        {
                            s += bl.ToString() + ",";
                        }
                        blcount++;
                    }
                    tw.WriteLine(s);
                }
                tw.Close();
            }
        }
    }
}
