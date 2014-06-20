using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Xml;
using System.IO;
using DOEgbXML;
using VectorMath;
using log4net;
using System.Text.RegularExpressions;

namespace gbXMLValidator
{
    class Program
    {
        //define tolerances for the tests
        static double coordtol = DOEgbXMLBasics.Tolerances.coordToleranceIP;
        //path of Text file output
        static string outputpath = @"C:\gbXML";
        static string reportpath = @"C:\gbXML\report.txt";
        //a good time to figure out what to do with my units

        //an example main that has all of the tests that we would like to perform
        static void Main(string[] args)
        {
            bool spaceBoundsPresent = false;
            //Destroy old reporting files
            if(File.Exists(reportpath))
            {
                File.Delete(reportpath);
            }
            //Basic File Preparation and Checks --------------------------------------------------------------------------------------------
            
            DOEgbXML.XMLParser parser = new XMLParser();
            //1-load the file
            //2-check for valid XML and valid gbXML against the XSD

            //hardcoded for testing purposes, eventually this file will be uploaded, or sent via a Restful API call
            //string path = @"C:\Users\Chiensi\Documents\C\CarmelSoft\gbXML Project\Phase 2\Validator Test Cases\Test Case 2.xml";
            //string path = @"C:\gbXML\test.xml";
            string path = @"C:\Users\Chiensi\Documents\C\CarmelSoft\gbXML Project\Phase 2\Validator Test Cases\Test Case 1 - Standard File - Three AdjacentSpaceIds.xml";
            XmlReader xmlreader = XmlReader.Create(path);
            XmlDocument myxml = new XmlDocument();
            myxml.Load(xmlreader);

            //3-get the namespace
            XmlNamespaceManager nsm = parser.getnsmanager(myxml);
            //figure out if metric or USIP (we have not found a reason to use this yet)
            parser.getunits(nsm, myxml);
            //Begin Parsing the XML and reporting on it----------------------------------------------------------
            //make a reporting object
            DOEgbXMLPhase2Report report = new DOEgbXMLPhase2Report();

            //Basic Uniqueness Constraints check-------------------------------------------------

            //Basic requirements check


            //ensure that all names of spaces are unique
            report.testType = TestType.Unique_Space_ID_Test;
            report = DOEgbXML.gbXMLSpaces.UniqueSpaceIdTest2(myxml, nsm, report);
            //process report
            ProcessReport(report,reportpath);
            report.Clear();
            
            XmlNodeList nodes = myxml.SelectNodes("/gbXMLv5:gbXML/gbXMLv5:Campus/gbXMLv5:Building/gbXMLv5:Space/gbXMLv5:SpaceBoundary", nsm);
            if (nodes.Count > 0)
            {
                spaceBoundsPresent = true;
                report.testType = TestType.Unique_Space_Boundary;
                //tests to ensure that each Space element only has the space boundary id called out once.
                report = DOEgbXML.gbXMLSpaces.UniqueSpaceBoundaryIdTest2(myxml, nsm, report);
                //process report
                ProcessReport(report, reportpath);
                report.Clear();
            }
            else
            {
                //needs to be included so the report can be processed
                report.testType = TestType.Unique_Space_Boundary;
                report.passOrFail = true;
                report.longMsg = "A test is usually performed here to ensure Space Boundaries have valid naming conventions.  This test was skipped (legally) because your file does not have space boundaries present.  Continuing to next test.";
                ProcessReport(report, reportpath);
            }

            //Space Tests
            //make a simplified representation of the spaces
            List<DOEgbXML.gbXMLSpaces> spaces = DOEgbXML.gbXMLSpaces.getSimpleSpaces(myxml, nsm);



            //4-check for self-intersecting polygons
            //report = DOEgbXML.gbXMLSpaces.SpaceSurfacesSelfIntersectionTest(spaces, report);
            //process report
            
            report.Clear();
            //check that all polyloops are in a counterclockwise direction

            report = DOEgbXML.gbXMLSpaces.SpaceSurfacesCCTest2(spaces, report);
            report.testType = TestType.Space_Surfaces_CC;
            //process report
            ProcessReport(report,reportpath);
            report.Clear();
            //-check for non-planar objects for all Spaces' polyloops
            report.testType = TestType.Space_Surfaces_Planar;
            report.coordtol = DOEgbXMLBasics.Tolerances.VectorAngleTolerance;
            report = DOEgbXML.gbXMLSpaces.SpaceSurfacesPlanarTest(spaces, report);
            //process report
            ProcessReport(report,reportpath);
            report.Clear();


            //valid space enclosure?
            report.tolerance = 0.0001;
            //when we are comparing angles in this function, we are testing the angle between dot products
            report.vectorangletol = DOEgbXMLBasics.Tolerances.dotproducttol;
            report.lengthtol = DOEgbXMLBasics.Tolerances.lengthTolerance;
            //toler
            report.coordtol = DOEgbXMLBasics.Tolerances.coordToleranceIP;
            report = CheckSpaceEnclosureSG(spaces, report);
            //process report
            ProcessReport(report, reportpath);
            report.Clear();


            //Surface tests----------------------------------------------------------------------------------
            //Basic Requirements ------------------------------------------------------

            //ensure that all names of surfaces are unique
            report.testType = TestType.Surface_ID_Uniqueness;
            report = DOEgbXML.SurfaceDefinitions.SurfaceIDUniquenessTest(myxml, nsm, report);
            //process report
            ProcessReport(report, reportpath);
            report.Clear();

            //Are there at least 4 surface definitions?  (see the surface requirements at the campus node)
            report.testType = TestType.At_Least_4_Surfaces;
            report = SurfaceDefinitions.AtLeast4Surfaces(myxml, nsm, report);
            //process report
            ProcessReport(report, reportpath);
            report.Clear();
            //Does the AdjacentSpaceId not exceed the max number allowable?
            //this needs to be updated!
            report.testType = TestType.Two_Adj_Space_Id;
            report = SurfaceDefinitions.AtMost2SpaceAdjId(myxml, nsm, report);
            //process report
            ProcessReport(report, reportpath);
            report.Clear();
            //Are all required elements and attributes in place?
            //report.testType = TestType.Required_Surface_Fields;
            //report = SurfaceDefinitions.RequiredSurfaceFields(myxml, nsm, report);
            //process report
            //ProcessReport(report, reportpath);
            //report.Clear();


            //now grab all the surfaceIds and make them available
            List<string> spaceIds = new List<string>();
            foreach(gbXMLSpaces s in spaces)
            {
                spaceIds.Add(s.id);
            }
            
            List<SurfaceDefinitions> surfaces = DOEgbXML.XMLParser.MakeSurfaceList(myxml, nsm);
            //make sure the surface Adjacent space Id names match only the the space Ids gathered above.  The adjacent space Ids can't have their own special values
            report.testType = TestType.Surface_Adj_Id_Match;
            report = DOEgbXML.SurfaceDefinitions.SurfaceAdjSpaceIdTest(spaceIds,surfaces, report);
            ProcessReport(report,reportpath);
            report.Clear();

            if (spaceBoundsPresent)
            {
                report.tolerance = DOEgbXMLBasics.Tolerances.coordToleranceIP;
                report.testType = TestType.Surface_ID_SB_Match;
                report = SurfaceDefinitions.SurfaceMatchesSpaceBoundary(myxml, nsm, report);
                //process report
                ProcessReport(report, reportpath);
                report.Clear();
            }
            else
            {
                report.testType = TestType.Surface_ID_SB_Match;
                report.passOrFail = true;
                report.longMsg = "A test is usually performed here to ensure Space Boundaries and Surfaces share the same ID.  This test was skipped (legally) because your file does not have space boundaries present.  Continuing to next test.";
                ProcessReport(report, reportpath);
            }

            //Does the polyloop right hand rule vector form the proper azimuth and tilt? (with and without a CADModelAzimuth)
            report.tolerance = DOEgbXMLBasics.Tolerances.VectorAngleTolerance;
            report.vectorangletol = DOEgbXMLBasics.Tolerances.VectorAngleTolerance;
            report.testType = TestType.Surface_Tilt_Az_Check;
            report = SurfaceDefinitions.SurfaceTiltAndAzCheck(myxml, nsm, report);
            //process report
            ProcessReport(report, reportpath);
            report.Clear();

            //planar surface test
            report.vectorangletol = DOEgbXMLBasics.Tolerances.VectorAngleTolerance;
            report = SurfaceDefinitions.TestSurfacePlanarTest(surfaces, report);
            //process report
            ProcessReport(report, reportpath);
            report.Clear();

            //I must take the surfaces, group them, and rearrange any interior surfaces' coordinates that should be pointed the opposite way
            string searchpath = "/gbXMLv5:gbXML/gbXMLv5:Campus/gbXMLv5:Building/gbXMLv5:Space";
            List<string> spaceids = DOEgbXML.gbXMLSpaces.getSpaceIds(myxml, nsm, searchpath);
            Dictionary<string, List<SurfaceDefinitions>> enclosure = new Dictionary<string, List<SurfaceDefinitions>>();
            foreach (string id in spaceids)
            {
                //find all surfaces with this adjacent space id
                //get their polyloops
                //and then match their polyloops
                List<SurfaceDefinitions> surflist = new List<SurfaceDefinitions>();
                foreach (SurfaceDefinitions surf in surfaces)
                {
                    
                    foreach (var adj in surf.AdjSpaceId)
                    {
                        if (adj == id)
                        {
                            surflist.Add(surf);
                            break;
                        }

                    }
                }
                enclosure[id] = surflist;
            }

            //counter clockwise winding test
            report.testType = TestType.Surface_CC_Test;
            report = SurfaceDefinitions.SurfaceCCTest(enclosure, report);
            //process report
            ProcessReport(report, reportpath);
            report.Clear();

            //self intersecting polygon test
            //report = SurfaceDefinitions.SurfaceSelfIntersectionTest(surfaces, report);
            //process the report
            report.Clear();


            //Is the Lower Left Corner properly defined?

            //surface enclosure tests
            report.tolerance = 0.0001;
            report.vectorangletol = DOEgbXMLBasics.Tolerances.VectorAngleTolerance;
            report.lengthtol = DOEgbXMLBasics.Tolerances.lengthTolerance;
            report.coordtol = DOEgbXMLBasics.Tolerances.coordToleranceIP;
            report.testType = TestType.Check_Surface_Enclosure;
            report = CheckSurfaceEnclosure(enclosure, report);
            ProcessReport(report, reportpath);
            report.Clear();

          
            //Openings Tests-----------------------------------------------------

            //Shading Devices Tests----------------------------------------------


            

        }

        private static void ProcessReport(DOEgbXMLPhase2Report report, string path)
        {
            if (!File.Exists(path))
            {


                using (System.IO.StreamWriter file = new System.IO.StreamWriter(path))
                {
                    file.WriteLine("Explanation of Test: " + report.testSummary);
                    if (report.passOrFail) { file.WriteLine("Test has Passed."); }
                    else { file.WriteLine("Test has Failed."); }
                    file.WriteLine("Explanation of Why Test Passed or Failed: " + report.longMsg);
                    if (report.TestPassedDict.Count() > 0)
                    {
                        file.WriteLine("Summary of findings: ");
                        foreach (KeyValuePair<string, bool> kp in report.TestPassedDict)
                        {
                            // If the line doesn't contain the word 'Second', write the line to the file. 
                            string line = kp.Key + ":" + kp.Value.ToString();
                            file.WriteLine(line);
                        }

                    }
                    //more detail


                    if (report.MessageList.Count() > 0)
                    {
                        foreach (KeyValuePair<string, List<string>> message in report.MessageList)
                        {
                            
                            foreach (string finding in message.Value)
                            {
                                string line = message.Key + ": ";
                                line += finding;
                                file.WriteLine(line);
                            }
                        }
                    }
                    file.WriteLine("\n");
                }
            }
            else
            {
                using (StreamWriter sw = File.AppendText(path))
                {
                    sw.WriteLine("Explanation of Test: " + report.testSummary);
                    if (report.passOrFail) { sw.WriteLine("Test has Passed."); }
                    else { sw.WriteLine("Test has failed."); }
                    sw.WriteLine("Explanation of Why Test Passed or Failed: " + report.longMsg);
                    if (report.TestPassedDict.Count() > 0)
                    {
                        sw.WriteLine("Summary of findings: ");
                        foreach (KeyValuePair<string, bool> kp in report.TestPassedDict)
                        {
                            // If the line doesn't contain the word 'Second', write the line to the file. 
                            string line = kp.Key + ":" + kp.Value.ToString();
                            sw.WriteLine(line);
                        }

                    }
                    //more detail

                    if (report.MessageList.Count() > 0)
                    {
                        foreach (KeyValuePair<string, List<string>> message in report.MessageList)
                        {
                            
                            foreach (string finding in message.Value)
                            {
                                string line = message.Key + ": ";
                                line += finding;
                                sw.WriteLine(line);
                            }
                        }
                    }
                    sw.WriteLine("\n");
                }
            }
        }
        //April 14, 2014 - Deprecated as not useful or within scope
        //private static Dictionary<Vector.CartCoord, Tuple<List<string>, List<bool>>> GetSBVertices(XmlNamespaceManager nsm, XmlDocument zexml)
        //{
        //    Dictionary<Vector.CartCoord,Tuple<List<string>,List<bool>>> sbdict = new Dictionary<Vector.CartCoord, Tuple<List<string>,List<bool>>>();

        //    //get each space boundary
        //    XmlNodeList nodes = zexml.SelectNodes("/gbXMLv5:gbXML/gbXMLv5:Campus/gbXMLv5:Building/gbXMLv5:Space/gbXMLv5:SpaceBoundary", nsm);

        //    int sbcount = 0;
        //    foreach (XmlNode sbnode in nodes)
        //    {
        //        //record the name of the space boundary
        //        string sbname="";
        //        XmlAttributeCollection spaceAtts = sbnode.Attributes;
        //        foreach (XmlAttribute at in spaceAtts)
        //        {
        //            if (at.Name == "surfaceIdRef")
        //            {
        //                sbname = at.Value;
        //                break;
        //            }
        //        }
        //        //create list of unique vertices
        //        XmlNode closedshell = sbnode.FirstChild;
        //        XmlNode PL = closedshell.FirstChild;

        //        foreach (XmlNode cp in PL)
        //        {
                    
        //            if (cp.Name == "CartesianPoint")
        //            {
        //                sbdict = AddCoordinate(sbcount, sbname, sbdict, cp);
        //            }
        //        }
        //        //XmlNodeList cnodes = sbnode.SelectNodes("/gbXMLv5:gbXML/gbXMLv5:Campus/gbXMLv5:Building/gbXMLv5:Space/gbXMLv5:SpaceBoundary/gbXMLv5:PlanarGeometry//gbXMLv5:PolyLoop", nsm);
        //        sbcount++;
        //    }
        //    return sbdict;
        //}

        //private static Dictionary<Vector.CartCoord, Tuple<List<string>, List<bool>>> GetSurfVertices(XmlNamespaceManager nsm, XmlDocument zexml)
        //{
        //    Dictionary<Vector.CartCoord, Tuple<List<string>, List<bool>>> surfdict = new Dictionary<Vector.CartCoord, Tuple<List<string>, List<bool>>>();

        //    //get each space boundary
        //    XmlNodeList nodes = zexml.SelectNodes("/gbXMLv5:gbXML/gbXMLv5:Campus/gbXMLv5:Surface", nsm);

        //    int surfcount = 0;
        //    foreach (XmlNode surfnode in nodes)
        //    {
        //        //record the name of the space boundary
        //        string surfname = "";
        //        XmlAttributeCollection spaceAtts = surfnode.Attributes;
        //        foreach (XmlAttribute at in spaceAtts)
        //        {
        //            if (at.Name == "id")
        //            {
        //                surfname = at.Value;
        //                break;
        //            }
        //        }
        //        //create list of unique vertices
        //        foreach (XmlNode childnode in surfnode)
        //        {
        //            if (childnode.Name == "PlanarGeometry")
        //            {
        //                XmlNode PL = childnode.FirstChild;
        //                foreach (XmlNode cp in PL)
        //                {

        //                    if (cp.Name == "CartesianPoint")
        //                    {
        //                        surfdict = AddCoordinate(surfcount, surfname, surfdict, cp);
        //                    }
        //                }
        //                surfcount++;
        //            }
        //        }
        //    }
        //    return surfdict;
        //}

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

        public static Dictionary<int, DOEgbXMLBasics.EdgeFamily> FindEdgeMatch(Dictionary<int, DOEgbXMLBasics.EdgeFamily> uniqueedges, DOEgbXMLBasics.EdgeFamily edge)
        {
            try
            {
                int edgeloopcounter = 0; //keeps track of how many guest edges in the unique edge dictionary I've searched through
                foreach (KeyValuePair<int, DOEgbXMLBasics.EdgeFamily> kp in uniqueedges)
                {
                    Vector.MemorySafe_CartCoord gueststartpt = kp.Value.startendpt[0];
                    //In the unique edge dictionary, I have located at least one point that is similar to my test edge, a tolerance needed?
                    if (gueststartpt.X == edge.startendpt[0].X && gueststartpt.Y == edge.startendpt[0].Y && gueststartpt.Z == edge.startendpt[0].Z)
                    {
                        //found at least one perfect coordinate match, try to match the second...get the endpoint of the unique edge in the dictionary
                        Vector.MemorySafe_CartCoord guestendpt = kp.Value.startendpt[1];
                        if (guestendpt.X == edge.startendpt[1].X && guestendpt.Y == edge.startendpt[1].Y && guestendpt.Z == edge.startendpt[1].Z)
                        {
                            //both match, means the match is perfect, so the unique edge has found it complement related edge.  Great!
                            kp.Value.relatedEdges.Add(edge);
                            //I am done searching this test edge, and I can start over again with the next edge in question
                            break;

                        }
                        else
                        {
                            //so far, I have found only one thing in common, sharing of one point, even though second point did not match, the edges could still align
                            //draw vector A
                            double Ax = guestendpt.X - edge.startendpt[1].X;
                            double Ay = guestendpt.Y - edge.startendpt[1].Y;
                            double Az = guestendpt.Z - edge.startendpt[1].Z;
                            Vector.MemorySafe_CartVect A = new Vector.MemorySafe_CartVect(Ax, Ay, Az);
                            double Amag = Vector.VectorMagnitude(A);

                            //take cross product to see if they are even in same plane
                            double evX = guestendpt.X - gueststartpt.X;
                            double evY = guestendpt.Y - gueststartpt.Y;
                            double evZ = guestendpt.Z - gueststartpt.Z;
                            Vector.MemorySafe_CartVect ev = new Vector.MemorySafe_CartVect(evX, evY, evZ);
                            double evmag = Vector.VectorMagnitude(ev);
                            Vector.MemorySafe_CartVect cross = Vector.CrossProduct(A, ev);
                            double dot = Vector.DotProduct(A, ev);
                            double crossmag = Vector.VectorMagnitude(cross);
                            //If Vector A and ev are parallel, then the cross product magnitude should be zero, add a small tolerance?
                            if (dot == 1)
                            {
                                //then we are at least parallel but they are perfect matches
                                //now see if both points of the test edge resides ON the  guest edge or OUTSIDE of it
                                double Bx = gueststartpt.X - edge.startendpt[1].X;
                                double By = gueststartpt.Y - edge.startendpt[1].Y;
                                double Bz = gueststartpt.Z - edge.startendpt[1].Z;
                                Vector.MemorySafe_CartVect B = new Vector.MemorySafe_CartVect(Bx, By, Bz);
                                double Bmag = Vector.VectorMagnitude(B);
                                //check to see if the test edge is inside the guest edge (shorter than)
                                //check to see if the test edge is outside of the guest edge (longer than)
                                //the first check is easier and should always yield an easy answer
                                if (Amag < evmag && Bmag < evmag)
                                {
                                    //true
                                    //we choose to create related edges relationship for test and guest edge because a perfect match wasn't found
                                    //so this test edge is added to the guest unique edge
                                    kp.Value.relatedEdges.Add(edge);
                                    //and the test edge itself accumulate its own relationships (it is seen as sort of unique)
                                    edge.relatedEdges.Add(kp.Value);
                                    //we will continue checking the unique guest edges to look for another match, since we have not found perfection
                                    edgeloopcounter++;
                                    continue;
                                }
                                //otherwise it does not fall inside the guest edge
                                //check now to see if the test edge falls outside of the guest edge
                                //it either overlaps it, or does not overlap it at all
                                //we know the two vectors are parallel and not antiparallel because of the dot product above
                                //now test the quality of their relationship
                            }
                            //some sort of tolerance here?
                            else if (dot == -1)
                            {
                                double testedgeX = edge.startendpt[1].X - edge.startendpt[0].X;
                                double testedgeY = edge.startendpt[1].Y - edge.startendpt[0].Y;
                                double testedgeZ = edge.startendpt[1].Z - edge.startendpt[0].Z;
                                Vector.MemorySafe_CartVect edgevec = new Vector.MemorySafe_CartVect(testedgeX, testedgeY, testedgeZ);
                                double testedgemag = Vector.VectorMagnitude(edgevec);

                                if (evmag < testedgemag)
                                {
                                    //this means the test edge is longer than the guest edge,
                                    //since we know they are parallel and share a common first point, we can conclude
                                    //they overlap but are not perfect matches.
                                    //so each accumulates a relationship
                                    kp.Value.relatedEdges.Add(edge);
                                    //the edge is still unique but accumulates a neighbor
                                    edge.relatedEdges.Add(kp.Value);
                                    edgeloopcounter++;
                                    continue;
                                }
                                //deprecated as not productive code on March 30, 2014

                                //double Cx = gueststartpt.X - edge.startendpt[1].X;
                                //double Cy = gueststartpt.Y - edge.startendpt[1].Y;
                                //double Cz = gueststartpt.Z - edge.startendpt[1].Z;
                                //Vector.MemorySafe_CartVect C = new Vector.MemorySafe_CartVect(Cx, Cy, Cz);
                                //double Cmag = Vector.VectorMagnitude(C);

                                //double Dx = guestendpt.X - edge.startendpt[1].X;
                                //double Dy = guestendpt.Y - edge.startendpt[1].Y;
                                //double Dz = guestendpt.Z - edge.startendpt[1].Z;
                                //Vector.MemorySafe_CartVect D = new Vector.MemorySafe_CartVect(Dx, Dy, Dz);
                                //double Dmag = Vector.VectorMagnitude(D);

                                //if (Dmag < testedgemag && Cmag < testedgemag)
                                //{
                                //    //this means the test edge is longer than the guest edge,
                                //    //since we know they are parallel and share a common first point, we can conclude
                                //    //they overlap but are not perfect matches.
                                //    //so each accumulates a relationship
                                //    kp.Value.relatedEdges.Add(edge);
                                //    //the edge is still unique but accumulates a neighbor
                                //    edge.relatedEdges.Add(kp.Value);
                                //    edgeloopcounter++;
                                //    continue;
                                //}

                            }//are the two lines parallel
                            else
                            {
                                //the two lines aren't parallel, so just move on
                                edgeloopcounter++;
                                continue;
                            }
                        }//else there is not a perfect match
                    } //if one coordinate has perfectly matched
                    //here the the start point of the guest edge matches the end point of the testedge, a switcheroo
                    else if (gueststartpt.X == edge.startendpt[1].X && gueststartpt.Y == edge.startendpt[1].Y && gueststartpt.Z == edge.startendpt[1].Z)
                    {
                        //found at least one perfect coordinate match, try to match the second
                        Vector.MemorySafe_CartCoord guestendpt = kp.Value.startendpt[1];
                        if (guestendpt.X == edge.startendpt[0].X && guestendpt.Y == edge.startendpt[0].Y && guestendpt.Z == edge.startendpt[0].Z)
                        {
                            //both match, means the match is perfect, so add it to the related surfaces list
                            kp.Value.relatedEdges.Add(edge);
                            break;

                        }
                        else
                        {
                            //the edge may be unique, though it could still have neighboring relationships
                            double Ax = guestendpt.X - edge.startendpt[0].X;
                            double Ay = guestendpt.Y - edge.startendpt[0].Y;
                            double Az = guestendpt.Z - edge.startendpt[0].Z;
                            Vector.MemorySafe_CartVect A = new Vector.MemorySafe_CartVect(Ax, Ay, Az);
                            double Amag = Vector.VectorMagnitude(A);

                            //take cross product to see if they are even in same plane
                            double evX = guestendpt.X - gueststartpt.X;
                            double evY = guestendpt.Y - gueststartpt.Y;
                            double evZ = guestendpt.Z - gueststartpt.Z;
                            Vector.MemorySafe_CartVect ev = new Vector.MemorySafe_CartVect(evX, evY, evZ);
                            double evmag = Vector.VectorMagnitude(ev);
                            double dot = Vector.DotProduct(A, ev);
                            Vector.MemorySafe_CartVect cross = Vector.CrossProduct(A, ev);
                            double crossmag = Vector.VectorMagnitude(cross);
                            //tolerance?
                            //we know that they are parallel, and therefore that the test edge is shorter than the guest edge
                            if (dot == 1)
                            {
                                //we now verify if the point resides on the edge or outside of it
                                double Bx = gueststartpt.X - edge.startendpt[0].X;
                                double By = gueststartpt.Y - edge.startendpt[0].Y;
                                double Bz = gueststartpt.Z - edge.startendpt[0].Z;
                                Vector.MemorySafe_CartVect B = new Vector.MemorySafe_CartVect(Bx, By, Bz);
                                double Bmag = Vector.VectorMagnitude(B);
                                //check to see if the test edge is inside the guest edge
                                if (Amag < evmag && Bmag < evmag)
                                {
                                    //this means it lies on the plane at least, so it shares, but it is also still independent because a perfect match wasn't found
                                    kp.Value.relatedEdges.Add(edge);
                                    //accumulate its own relationships
                                    edge.relatedEdges.Add(kp.Value);
                                    edgeloopcounter++;
                                    continue;
                                }
                            }
                            //we know the lines are antiparallel, meaning that the test edge is longer than the guest edge
                            else if (dot == -1)
                            {
                                double testedgeX = edge.startendpt[1].X - edge.startendpt[0].X;
                                double testedgeY = edge.startendpt[1].Y - edge.startendpt[0].Y;
                                double testedgeZ = edge.startendpt[1].Z - edge.startendpt[0].Z;
                                Vector.MemorySafe_CartVect edgevec = new Vector.MemorySafe_CartVect(testedgeX, testedgeY, testedgeZ);
                                double testedgemag = Vector.VectorMagnitude(edgevec);

                                //we verify
                                if (testedgemag > evmag)
                                {
                                    //this means the test edge is longer than the guest edge, but they overlap
                                    kp.Value.relatedEdges.Add(edge);
                                    //the edge is still unique but accumulates a neighbor
                                    edge.relatedEdges.Add(kp.Value);
                                    edgeloopcounter++;
                                    continue;
                                }
                                //deprecated as not productive code March 30 2014
                                //double Cx = gueststartpt.X - edge.startendpt[0].X;
                                //double Cy = gueststartpt.Y - edge.startendpt[0].Y;
                                //double Cz = gueststartpt.Z - edge.startendpt[0].Z;
                                //Vector.MemorySafe_CartVect C = new Vector.MemorySafe_CartVect(Cx, Cy, Cz);
                                //double Cmag = Vector.VectorMagnitude(C);

                                //double Dx = guestendpt.X - edge.startendpt[0].X;
                                //double Dy = guestendpt.Y - edge.startendpt[0].Y;
                                //double Dz = guestendpt.Z - edge.startendpt[0].Z;
                                //Vector.MemorySafe_CartVect D = new Vector.MemorySafe_CartVect(Dx, Dy, Dz);
                                //double Dmag = Vector.VectorMagnitude(D);

                                //if (Dmag < testedgemag && Cmag < testedgemag)
                                //{
                                //    //this means the test edge is longer than the guest edge, but they overlap
                                //    kp.Value.relatedEdges.Add(edge);
                                //    //the edge is still unique but accumulates a neighbor
                                //    edge.relatedEdges.Add(kp.Value);
                                //    edgeloopcounter++;
                                //    continue;
                                //}
                            }
                            else
                            {
                                //this other point isn't relevant, and the edges don't coincide
                                edgeloopcounter++;
                                continue;
                            }
                        }

                    }//second point matches first point
                    //neither points perfectly coincide, so we do an exhaustive overlap check.
                    else
                    {
                        Vector.MemorySafe_CartCoord guestendpt = kp.Value.startendpt[1];
                        //are the two vectors even parallel?  because if they are not, no need to get more complex
                        double evX = guestendpt.X - gueststartpt.X;
                        double evY = guestendpt.Y - gueststartpt.Y;
                        double evZ = guestendpt.Z - gueststartpt.Z;
                        Vector.MemorySafe_CartVect guestvec = new Vector.MemorySafe_CartVect(evX, evY, evZ);
                        double guestvecmag = Vector.VectorMagnitude(guestvec);

                        double edgeX = edge.startendpt[1].X - edge.startendpt[0].X;
                        double edgeY = edge.startendpt[1].Y - edge.startendpt[0].Y;
                        double edgeZ = edge.startendpt[1].Z - edge.startendpt[0].Z;
                        Vector.MemorySafe_CartVect testedgev = new Vector.MemorySafe_CartVect(edgeX, edgeY, edgeZ);
                        //tolerance?
                        if (Vector.VectorMagnitude(Vector.CrossProduct(guestvec, testedgev)) != 0)
                        {
                            //they are not even parallel so move on
                            edgeloopcounter++;
                            continue;
                        }
                        //if the cross product is zero
                        //try to determine if the two edges are parallel or antiparallel
                        //test edge point 1

                        double Ax = gueststartpt.X - edge.startendpt[0].X;
                        double Ay = gueststartpt.Y - edge.startendpt[0].Y;
                        double Az = gueststartpt.Z - edge.startendpt[0].Z;
                        Vector.MemorySafe_CartVect A = new Vector.MemorySafe_CartVect(Ax, Ay, Az);
                        double Amag = Vector.VectorMagnitude(A);

                        //take cross product to see if they are even in same plane

                        Vector.MemorySafe_CartVect cross1 = Vector.CrossProduct(A, guestvec);
                        double crossA = Vector.VectorMagnitude(cross1);
                        //tolerance?
                        if (crossA == 0)
                        {
                            //we are at least parallel, now to check for a real intersection
                            double Bx = gueststartpt.X - edge.startendpt[1].X;
                            double By = gueststartpt.Y - edge.startendpt[1].Y;
                            double Bz = gueststartpt.Z - edge.startendpt[1].Z;
                            Vector.MemorySafe_CartVect B = new Vector.MemorySafe_CartVect(Bx, By, Bz);
                            double Bmag = Vector.VectorMagnitude(B);
                            double crossB = Vector.VectorMagnitude(Vector.CrossProduct(B, guestvec));
                            //check to see if the test edge's first point (index 0) is totally inside the guest edge
                            if (crossB == 0)
                            {
                                //we know they are now surely parallel, and that they are on top of one another, question is, do they overlap?
                                Vector.PointVector guestPV = new Vector.PointVector(gueststartpt, guestendpt);
                                Vector.PointVector test1 = new Vector.PointVector(gueststartpt, edge.startendpt[0]);
                                Vector.PointVector test2 = new Vector.PointVector(gueststartpt, edge.startendpt[1]);

                                bool test1inside = false;
                                bool test2inside = false;
                                //tolerance?
                                if (Vector.DotProduct(guestPV.v, test1.v) == 1)
                                {
                                    //then pointed in the same direction

                                    double ratio = Amag / guestvecmag;
                                    if (ratio > 1)
                                    {
                                        //then this point is outside of evmag
                                    }
                                    else
                                    {
                                        //there is an intersection at the point of the test edge
                                        test1inside = true;
                                    }
                                }
                                //tolerance?
                                if (Vector.DotProduct(guestPV.v, test1.v) == 1)
                                {
                                    //then this second test edge point is in a straight line in the same direction as the vector formed by the test edge
                                    double ratio = Bmag / guestvecmag;
                                    if (ratio > 1)
                                    {
                                        //then this point is outside
                                    }
                                    else
                                    {
                                        //this is also inside
                                        test2inside = true;
                                    }

                                }

                                if (test1inside == true && test2inside == false)
                                {
                                    //make a temporary edge from just point 1 in the test edge
                                    //this should never happen for a single space, should it?

                                }
                                else if (test1inside == false && test2inside == true)
                                {
                                    //make a temporary edge from just point 2 in the test edge
                                    //this should never happen for a single space, should it?
                                }
                                else if (test1inside == true && test2inside == true)
                                {
                                    //the entire test edge is contained within
                                    //this means the test edge is longer than the guest edge, but they overlap
                                    kp.Value.relatedEdges.Add(edge);
                                    //the edge is still unique but accumulates a neighbor
                                    edge.relatedEdges.Add(kp.Value);
                                    edgeloopcounter++;
                                    continue;
                                }
                                else
                                {
                                    //do nothing because neither points reside inside.
                                }

                            }
                        }
                    }
                    int uniqueedgect = uniqueedges.Count() + 1;
                    uniqueedges[uniqueedgect] = edge;
                }

            }
            catch (Exception e)
            {

            }
            return uniqueedges;
        }

        public static DOEgbXMLPhase2Report CheckSpaceEnclosureSB(List<gbXMLSpaces> spaces, DOEgbXMLPhase2Report report)
        {
            try
            {
                report.passOrFail = true;
                foreach (gbXMLSpaces space in spaces)
                {
                    List<string> ml = new List<string>();
                    if (space.spacebounds.Count() > 0)
                    {
                        
                        ml.Add("Testing begins for Space Boundary water tightness.");
                        Dictionary<int, Vector.EdgeFamily> uniqueedges = new Dictionary<int, Vector.EdgeFamily>();
                        ml.Add(space.id + ": has SpaceBoundary representation.");
                        ml.Add(space.id + ": START checking space boundary enclosure.");
                        foreach (DOEgbXML.gbXMLSpaces.SpaceBoundary sb in space.spacebounds)
                        {
                            uniqueedges = Vector.GetEdgeFamilies(sb.surfaceIdRef, uniqueedges, sb.sbplane.pl.plcoords, .0001, .0001);
                        }
                        ml.Add("Gathered space boundary edges and neighboring relationships.");
                        //see how well enclosure is formed
                        ml.Add("Validating space boundary edge alignment with one another - water tightness check.");

                        //new function added April 11, 2014
                        report = MatchEdges(uniqueedges, ml, report, space.id);
                           
                    }
                    else
                    {
                        ml.Add(space.id + ": Has no Valid Space Boundaries.  This is not an error.  The validator will continue searching for other enclosed boundaries.");
                    }
                    report.MessageList[space.id] = ml;
                }
            }
            catch (Exception e)
            {
                report.longMsg = ("SORRY, we have run into an unexpected issue:" + e.ToString());
                report.passOrFail = false;
                return report;
            }
            return report;
        }

        public static DOEgbXMLPhase2Report CheckSpaceEnclosureSG(List<gbXMLSpaces> spaces, DOEgbXMLPhase2Report report)
        {
            try
            {
                report.testSummary = "This test checks the enclosure defined by the ShellGeometry PolyLoops for each given space in your gbXML file.";
                report.testSummary += " This is an optional test because ShellGeometry definitions are optional.";
                report.testSummary += " An enclosure test is important because it ensures that each of the surfaces in the gbXML definition is properly aligned ";
                report.testSummary += " with its neighbor.  The test checks to make sure that all edges of each surface line up with one another so that there are not";
                report.testSummary += "any gaps.";
                report.passOrFail = true;
                foreach (gbXMLSpaces space in spaces)
                {
                    List<string> ml = new List<string>();
                    Dictionary<int, Vector.EdgeFamily> uniqueedges = new Dictionary<int, Vector.EdgeFamily>();

                    if (space.sg.cs.ploops.Count() > 0)
                    {
                        string rep = space.id + " has ShellGeometry PolyLoops.  Conducting tests of ShellGeometry PolyLoops water tightness";
                        report.TestPassedDict[rep] = true;
                        int sgcount = 1;
                        foreach (DOEgbXML.gbXMLSpaces.PolyLoop pl in space.sg.cs.ploops)
                        {
                            string surfaceid = "shellgeometry-" + sgcount;
                            uniqueedges = Vector.GetEdgeFamilies(surfaceid, uniqueedges, pl.plcoords, report.coordtol, report.vectorangletol);
                            sgcount++;
                        }
                        if (uniqueedges.Count > 0)
                        {
                            string erep = space.id + ": Gathered edges of the ShellGeometry PolyLoop successfullly.";
                            report.TestPassedDict[erep] = true;
                        }
                        else
                        {
                            string erep = space.id + ": Gathered edges of the ShellGeometry PolyLoop unsuccessfullly.";
                            report.TestPassedDict[erep] = false;
                        }
                        
                        //see how well enclosure is formed
                        //new function added April 11, 2014
                        report = MatchEdges(uniqueedges, ml, report, space.id);
                    }
                    else
                    {
                        string rep = space.id + " does not has ShellGeometry PolyLoops (this is not an error).  Conducting tests of ShellGeometry PolyLoops water tightness";
                        report.TestPassedDict[rep] = false;
                    }
                    report.MessageList[space.id] = ml;
                }

            }
            catch (Exception e)
            {
                report.longMsg = ("SORRY, we have run into an unexpected issue:" + e.ToString());
                report.passOrFail = false;
                return report;
            }
            return report;
        }

        //April 11, 2014
        //This is the algorithm that attempts to find matches of edges on the enclosure.  It is used by all the Check Enclosure routines.
        //by Chien Harriman - Carmel Software Corporation
        public static DOEgbXMLPhase2Report MatchEdges(Dictionary<int,Vector.EdgeFamily> uniqueedges, List<string> ml, DOEgbXMLPhase2Report report, string spaceid)
        {
            try
            {
                int totaledgect = 0;
                int matchededges = 0;
                string lastedgenm = "";
                int surfedgect = 0;
                foreach (KeyValuePair<int, Vector.EdgeFamily> edgekp in uniqueedges)
                {
                    //a way to count edges
                    if (edgekp.Value.sbdec != lastedgenm)
                    {
                        //reset
                        lastedgenm = edgekp.Value.sbdec;
                        surfedgect = 0;
                    }
                    //here is the easiest case where there is only one related edge
                    //we know this must be a perfect match, or entirely envelopes the edge 
                    if (edgekp.Value.relatedEdges.Count() == 1)
                    {
                        Vector.MemorySafe_CartCoord edgestart = edgekp.Value.startendpt[0];
                        Vector.MemorySafe_CartCoord edgeend = edgekp.Value.startendpt[1];
                        Vector.MemorySafe_CartCoord relstart = edgekp.Value.relatedEdges[0].startendpt[0];
                        Vector.MemorySafe_CartCoord relend = edgekp.Value.relatedEdges[0].startendpt[1];
                        //if the lengths are the same, then they should match perfectly.
                        //this is a valid conclusion because we already have identified that they aligh and 
                        //are in the same space.
                        double edgeX = edgestart.X - edgeend.X;
                        double edgeY = edgestart.Y - edgeend.Y;
                        double edgeZ = edgestart.Z - edgeend.Z;
                        Vector.MemorySafe_CartVect edgev = new Vector.MemorySafe_CartVect(edgeX, edgeY, edgeZ);
                        double edgemag = Vector.VectorMagnitude(edgev);

                        double relX = relstart.X - relend.X;
                        double relY = relstart.Y - relend.Y;
                        double relZ = relstart.Z - relend.Z;
                        Vector.MemorySafe_CartVect relv = new Vector.MemorySafe_CartVect(relX, relY, relZ);
                        double relmag = Vector.VectorMagnitude(relv);
                        //do the check here to see if the two edges (current and related) are the same length
                        if (Math.Abs(relmag - edgemag) < report.coordtol)
                        {
                            
                            //should match perfectly
                            ml.Add(edgekp.Value.sbdec + " Edge " + surfedgect.ToString()+ " should have perfectly matched coordinates.");
                            List<bool> match = new List<bool>();
                            double tol = report.tolerance;
                            for (int i = 0; i < 2; i++)
                            {
                                Vector.MemorySafe_CartCoord p1 = edgekp.Value.relatedEdges[0].startendpt[i];
                                for (int j = 0; j < 2; j++)
                                {
                                    string x = p1.X.ToString();
                                    string y = p1.Y.ToString();
                                    string z = p1.Z.ToString();
                                    string coordstr = "(" + x + "," + y + "," + z + ")";
                                    Vector.MemorySafe_CartCoord p2 = edgekp.Value.startendpt[j];
                                    if (p2.X == p1.X && p2.Y == p1.Y && p2.Z == p1.Z)
                                    {
                                        match.Add(true);
                                        ml.Add("PERFECT MATCH: " +edgekp.Value.sbdec+ " Edge "+surfedgect.ToString()+" Coordinate "+ coordstr);
                                    }
                                    else if (Math.Abs(p2.X - p1.X) < tol && Math.Abs(p2.Y - p1.Y) < report.coordtol && Math.Abs(p2.Z - p1.Z) < report.coordtol)
                                    {
                                        match.Add(true);
                                        ml.Add("MATCH: " + edgekp.Value.sbdec + " Edge " + surfedgect.ToString() + " Coordinate " + coordstr);
                                        report.passOrFail = false;
                                    }
                                }
                            }
                            if (match.Count() == 2)
                            {
                                ml.Add("PASS: " + edgekp.Value.sbdec + " Edge " + surfedgect.ToString() + " forms a tight enclosure with its neighbor.");
                                report.passOrFail = true;
                                //we +2 here because the related edge is not recorded
                                totaledgect += 2;
                                matchededges += 2;
                            }
                            else
                            {
                                ml.Add("FAIL: " + edgekp.Value.sbdec + " Edge " + surfedgect.ToString() + " does not form a tight enclosure with its neighbor.");
                                report.passOrFail = false;
                                //we +2 here because the related edge is not recorded
                                totaledgect += 2;
                                matchededges += 0;
                            }
                        }
                        //April 7 2014
                        //it is safe to conclude that the the two are related, but they overlap.  In this case, since there is only one neighbor
                        //it should be the case that one edge entirely envelops the other edge.  
                        //this edge, has to be the one that is enveloped, because it only has one related edge, by convention
                        else
                        {
                            ml.Add(edgekp.Value.sbdec + " Edge " + surfedgect.ToString()+" should be enclosed by its neighboring edge.");
                            //es--------------ee
                            //rs-------------------re
                            if (Math.Abs(edgestart.X - relstart.X) <= report.coordtol && Math.Abs(edgestart.Y - relstart.Y) <= report.coordtol && Math.Abs(edgestart.Z - relstart.Z) <= coordtol)
                            {
                                string x = edgestart.X.ToString();
                                string y = edgestart.Y.ToString();
                                string z = edgestart.Z.ToString();
                                string coordstr = "(" + x + "," + y + "," + z + ")";
                                ml.Add("MATCH: " + edgekp.Value.sbdec + " Edge " + surfedgect.ToString() + " Coordinate " + coordstr);
                                double Cx = edgeend.X - relstart.X;
                                double Cy = edgeend.Y - relstart.Y;
                                double Cz = edgeend.Z - relstart.Z;
                                Vector.MemorySafe_CartVect C = new Vector.MemorySafe_CartVect(Cx, Cy, Cz);

                                double Dx = edgeend.X - relend.X;
                                double Dy = edgeend.Y - relend.Y;
                                double Dz = edgeend.Z - relend.Z;
                                Vector.MemorySafe_CartVect D = new Vector.MemorySafe_CartVect(Dx, Dy, Dz);

                                double dotend = Vector.DotProductMag(C, D);
                                //both of these dot products should point in opposite directions, proving the edge is entirely enveloped
                                if (Math.Abs(dotend) - 1 <= report.vectorangletol)
                                {
                                    ml.Add("PASS" + edgekp.Value.sbdec + " Edge " + surfedgect.ToString() + " overlapping neighbor forms a tight enclosure with its neighbor.");
                                }
                                else
                                {
                                    ml.Add("FAIL" + edgekp.Value.sbdec + " Edge " + surfedgect.ToString() + "overlapping neighbor does not form a tight enclosure with its neighbor.");
                                    report.passOrFail = false;
                                }
                            }
                            //es-----------------ee
                            //re--------------------------rs
                            else if (Math.Abs(edgestart.X - relend.X) <= report.coordtol && Math.Abs(edgestart.Y - relend.Y) <= report.coordtol && Math.Abs(edgestart.Z - relend.Z) <= coordtol)
                            {
                                string x = edgestart.X.ToString();
                                string y = edgestart.Y.ToString();
                                string z = edgestart.Z.ToString();
                                string coordstr = "(" + x + "," + y + "," + z + ")";
                                ml.Add("MATCH: " + edgekp.Value.sbdec + " Edge " + surfedgect.ToString() + " Coordinate " + coordstr);
                                double Cx = edgeend.X - relstart.X;
                                double Cy = edgeend.Y - relstart.Y;
                                double Cz = edgeend.Z - relstart.Z;
                                Vector.MemorySafe_CartVect C = new Vector.MemorySafe_CartVect(Cx, Cy, Cz);

                                double Dx = edgeend.X - relend.X;
                                double Dy = edgeend.Y - relend.Y;
                                double Dz = edgeend.Z - relend.Z;
                                Vector.MemorySafe_CartVect D = new Vector.MemorySafe_CartVect(Dx, Dy, Dz);

                                double dotend = Vector.DotProductMag(C, D);
                                //both of these dot products should point in opposite directions, proving the edge is entirely enveloped
                                if (Math.Abs(dotend) - 1 <= report.vectorangletol)
                                {
                                    ml.Add("PASS" + edgekp.Value.sbdec + " Edge " + surfedgect.ToString() + " overlapping neighbor forms a tight enclosure with its neighbor.");
                                    totaledgect += 1;
                                    matchededges += 1;
                                }
                                else
                                {
                                    ml.Add("FAIL" + edgekp.Value.sbdec + " Edge " + surfedgect.ToString() + " overlapping neighbor does not form a tight enclosure with its neighbor.");
                                    report.passOrFail = false;
                                    totaledgect += 1;
                                    matchededges += 0;
                                }
                            }
                            //ee-----------------es
                            //rs--------------------------re
                            else if (Math.Abs(edgeend.X - relstart.X) <= report.coordtol && Math.Abs(edgeend.Y - relstart.Y) <= report.coordtol && Math.Abs(edgeend.Z - relstart.Z) <= coordtol)
                            {
                                string x = edgeend.X.ToString();
                                string y = edgeend.Y.ToString();
                                string z = edgeend.Z.ToString();
                                string coordstr = "(" + x + "," + y + "," + z + ")";
                                ml.Add("MATCH: " + edgekp.Value.sbdec + " Edge " + surfedgect.ToString() + " Coordinate " + coordstr);
                                double Ax = edgestart.X - relstart.X;
                                double Ay = edgestart.Y - relstart.Y;
                                double Az = edgestart.Z - relstart.Z;
                                Vector.MemorySafe_CartVect A = new Vector.MemorySafe_CartVect(Ax, Ay, Az);

                                double Bx = edgestart.X - relend.X;
                                double By = edgestart.Y - relend.Y;
                                double Bz = edgestart.Z - relend.Z;
                                Vector.MemorySafe_CartVect B = new Vector.MemorySafe_CartVect(Bx, By, Bz);

                                double dotstart = Vector.DotProductMag(A, B);
                                //both of these dot products should point in opposite directions, proving the edge is entirely enveloped
                                if (Math.Abs(dotstart) - 1 <= report.vectorangletol)
                                {
                                    ml.Add("PASS" + edgekp.Value.sbdec + " Edge " + surfedgect.ToString() + " overlapping neighbor forms a tight enclosure with its neighbor.");
                                    totaledgect += 1;
                                    matchededges += 1;
                                }
                                else
                                {
                                    ml.Add("FAIL" + edgekp.Value.sbdec + " Edge " + surfedgect.ToString() + " overlapping neighbor does not form a tight enclosure with its neighbor.");
                                    report.passOrFail = false;
                                    totaledgect += 1;
                                    matchededges += 0;
                                }
                            }
                            //ee-----------------es
                            //re--------------------------rs
                            else if (Math.Abs(edgeend.X - relend.X) <= report.coordtol && Math.Abs(edgeend.Y - relend.Y) <= report.coordtol && Math.Abs(edgeend.Z - relend.Z) <= coordtol)
                            {
                                string x = edgeend.X.ToString();
                                string y = edgeend.Y.ToString();
                                string z = edgeend.Z.ToString();
                                string coordstr = "(" + x + "," + y + "," + z + ")";
                                ml.Add("MATCH: " + edgekp.Value.sbdec + " Edge " + surfedgect.ToString() + " Coordinate " + coordstr);
                                double Ax = edgestart.X - relstart.X;
                                double Ay = edgestart.Y - relstart.Y;
                                double Az = edgestart.Z - relstart.Z;
                                Vector.MemorySafe_CartVect A = new Vector.MemorySafe_CartVect(Ax, Ay, Az);

                                double Bx = edgestart.X - relend.X;
                                double By = edgestart.Y - relend.Y;
                                double Bz = edgestart.Z - relend.Z;
                                Vector.MemorySafe_CartVect B = new Vector.MemorySafe_CartVect(Bx, By, Bz);

                                double dotstart = Vector.DotProductMag(A, B);
                                //both of these dot products should point in opposite directions, proving the edge is entirely enveloped
                                if (Math.Abs(dotstart) - 1 <= report.vectorangletol)
                                {
                                    ml.Add("PASS" + edgekp.Value.sbdec + " Edge " + surfedgect.ToString() + " overlapping neighbor forms a tight enclosure with its neighbor.");
                                    totaledgect += 1;
                                    matchededges += 1;
                                }
                                else
                                {
                                    ml.Add("FAIL" + edgekp.Value.sbdec + " Edge " + surfedgect.ToString() + " overlapping neighbor does not form a tight enclosure with its neighbor.");
                                    report.passOrFail = false;
                                    totaledgect += 1;
                                    matchededges += 0;
                                }
                            }
                            //overlapping edges
                            else
                            {
                                ml.Add(edgekp.Value.sbdec + " Edge " + surfedgect.ToString() + " is overlapped at both ends by its neighboring edge.");
                                double Ax = edgestart.X - relstart.X;
                                double Ay = edgestart.Y - relstart.Y;
                                double Az = edgestart.Z - relstart.Z;
                                Vector.MemorySafe_CartVect A = new Vector.MemorySafe_CartVect(Ax, Ay, Az);

                                double Bx = edgestart.X - relend.X;
                                double By = edgestart.Y - relend.Y;
                                double Bz = edgestart.Z - relend.Z;
                                Vector.MemorySafe_CartVect B = new Vector.MemorySafe_CartVect(Bx, By, Bz);

                                double Cx = edgeend.X - relstart.X;
                                double Cy = edgeend.Y - relstart.Y;
                                double Cz = edgeend.Z - relstart.Z;
                                Vector.MemorySafe_CartVect C = new Vector.MemorySafe_CartVect(Cx, Cy, Cz);

                                double Dx = edgeend.X - relend.X;
                                double Dy = edgeend.Y - relend.Y;
                                double Dz = edgeend.Z - relend.Z;
                                Vector.MemorySafe_CartVect D = new Vector.MemorySafe_CartVect(Dx, Dy, Dz);

                                double dotstart = Vector.DotProductMag(A, B);
                                double dotend = Vector.DotProductMag(C, D);
                                //both of these dot products should point in opposite directions, proving the edge is entirely enveloped
                                if (Math.Abs(dotstart) - 1 <= report.vectorangletol && Math.Abs(dotend) - 1 <= report.vectorangletol)
                                {
                                    ml.Add("PASS" + edgekp.Value.sbdec + " Edge " + surfedgect.ToString() + " overlapping neighbor forms a tight enclosure with its neighbor.");
                                    totaledgect += 1;
                                    matchededges += 1;
                                }
                                else
                                {
                                    ml.Add("FAIL" + edgekp.Value.sbdec + " Edge " + surfedgect.ToString() + " overlapping neighbor does not form a tight enclosure with its neighbor.");
                                    report.passOrFail = false;
                                    totaledgect += 1;
                                    matchededges += 0;
                                }
                            }
                        }

                    }
                    else if (edgekp.Value.relatedEdges.Count() > 1)
                    {
                        ml.Add(edgekp.Value.sbdec + " Edge " + surfedgect.ToString() + " has " + edgekp.Value.relatedEdges.Count() + " neighboring edges.");
                        //more robust testing
                        Vector.EdgeFamily[] orderededges = new Vector.EdgeFamily[edgekp.Value.relatedEdges.Count()];
                        //align the related edges
                        orderededges = AlignEdges(orderededges, edgekp.Value, coordtol);
                        
                        int coordcount = 0;
                        double segmentslength = 0;
                        int lastct = edgekp.Value.relatedEdges.Count() * 2 - 2;
                        Vector.CartCoord st = new Vector.CartCoord();
                        Vector.CartCoord end = new Vector.CartCoord();
                        foreach (Vector.EdgeFamily edge in orderededges)
                        {

                            if (coordcount == 0)
                            {
                                st.X = edge.startendpt[0].X;
                                st.Y = edge.startendpt[0].Y;
                                st.Z = edge.startendpt[0].Z;
                            }

                            if (coordcount == lastct)
                            {
                                end.X = edge.startendpt[1].X;
                                end.Y = edge.startendpt[1].Y;
                                end.Z = edge.startendpt[1].Z;
                            }
                            Vector.MemorySafe_CartVect v = Vector.CreateMemorySafe_Vector(edge.startendpt[0], edge.startendpt[1]);
                            double mag = Vector.VectorMagnitude(v);
                            segmentslength += mag;
                            coordcount += 2;
                        }
                        Vector.MemorySafe_CartVect v2 = Vector.CreateMemorySafe_Vector(edgekp.Value.startendpt[0], edgekp.Value.startendpt[1]);
                        double mag2 = Vector.VectorMagnitude(v2);
                        if (Math.Abs(segmentslength - mag2) < report.lengthtol)
                        {
                            ml.Add(edgekp.Value.sbdec + ":PASS:Multiple Overlapping edges properly match.");
                        }
                        else
                        {
                            //then something is wrong.
                            ml.Add(edgekp.Value.sbdec + ":FAIL:Overlapping edges do not match as expected.");
                            report.passOrFail = false;
                        }

                    }
                    else if (edgekp.Value.relatedEdges.Count() == 0)
                    {
                        //something is wrong
                        ml.Add(edgekp.Value.sbdec+"FAIL:" + edgekp.Value.sbdec + " Edge " + surfedgect.ToString() + " has no reported neighboring edges.");
                        report.passOrFail = false;
                    }
                    surfedgect += 1;
                }
                report.MessageList[spaceid] = ml;
                //all edges align = true
                if (report.passOrFail)
                {
                    report.longMsg = "TEST PASSED: " + totaledgect + " edges in the gbXML file.  " + matchededges + " edges found with ideal alignment.";
                }
                //all edges align = false
                else
                {
                    report.longMsg = "TEST FAILED: " + totaledgect + " edges in the gbXML file.  " + matchededges + " edges found with ideal alignment.";
                }
                return report;

            }
            catch (Exception e)
            {
                report.longMsg = ("ERROR, we have run into an unexpected issue:" + e.ToString());
                report.passOrFail = false;
                return report;
            }

        }
        //we have already proven that the neighboring edges are aligned with the edge
        //so we do not perform an exhaustive vector search to find the answer.  All we do instead is try and position them correctly 
        //in sequence as well as is possible.  
        //we do not check for overlaps until later
        public static Vector.EdgeFamily[] AlignEdges(Vector.EdgeFamily[] alignededges, Vector.EdgeFamily edge, double tol)
        {
            try
            {
                List<double> magnitudes = new List<double>();
                List<int> indices = new List<int>();
                for (int i = 0; i < edge.relatedEdges.Count(); i++)
                {
                    Vector.MemorySafe_CartCoord longstart = edge.startendpt[0];
                    Vector.MemorySafe_CartCoord longend = edge.startendpt[1];
                    Vector.MemorySafe_CartVect longv = Vector.CreateMemorySafe_Vector(longstart, longend);
                    double maglong = Vector.VectorMagnitude(longv);
                    Vector.EdgeFamily shortedge = edge.relatedEdges[i];
                    Vector.MemorySafe_CartCoord shortstart = shortedge.startendpt[0];
                    Vector.MemorySafe_CartCoord shortend = shortedge.startendpt[1];
                    Vector.MemorySafe_CartVect el1 = Vector.CreateMemorySafe_Vector(longstart, shortstart);
                    Vector.MemorySafe_CartVect el2 = Vector.CreateMemorySafe_Vector(longstart, shortend);
                    double magel1 = Vector.VectorMagnitude(el1);
                    double magel2 = Vector.VectorMagnitude(el2);

                    
                    //put the greater of the two magnitudes in the list
                    if (magel1 > magel2)
                    {
                        shortedge.startendpt.Reverse();
                        if (magnitudes.Count >= 1)
                        {
                            bool added = false;
                            for (int m = 0; m < magnitudes.Count(); m++)
                            {
                                if (magel1 < magnitudes[m])
                                {
                                    magnitudes.Insert(m, magel1);
                                    indices.Insert(m, i);
                                    added = true;
                                    break;
                                }
                            }
                            if (!added)
                            {
                                magnitudes.Add(magel1);
                                indices.Add(i);
                            }

                        }
                        else
                        {
                            magnitudes.Add(magel1);
                            indices.Add(i);
                        }
                    }
                    else
                    {
                        if (magnitudes.Count >= 1)
                        {
                            bool added = false;
                            for (int m = 0; m < magnitudes.Count(); m++)
                            {
                                if (magel1 < magnitudes[m])
                                {
                                    magnitudes.Insert(m, magel2);
                                    indices.Insert(m, i);
                                    added = true;
                                    break;
                                }
                            }
                            if (!added)
                            {
                                magnitudes.Add(magel2);
                                indices.Add(i);
                            }

                        }
                        else
                        {
                            magnitudes.Add(magel2);
                            indices.Add(i);
                        }
                    }

                }
                int alignedcounter = 0;
                foreach (int i in indices)
                {
                    alignededges[alignedcounter] = edge.relatedEdges[i];
                    alignedcounter++;
                }

                return alignededges;

                
            }
            catch (Exception e)
            {

            }
            return alignededges;
        }

        public static DOEgbXMLPhase2Report CheckSurfaceEnclosure(Dictionary<string, List<SurfaceDefinitions>> surfaceEnclosures, DOEgbXMLPhase2Report report)
        {
            try
            {
                report.testSummary = "This test checks surfaces proclaiming to be children of a given space ID.";
                report.testSummary += "  It searches each of the surfaces' edges and tries to find other edges that align.";
                report.passOrFail = true;
                foreach (KeyValuePair<string, List<SurfaceDefinitions>> kp in surfaceEnclosures)
                {
                    List<string> ml = new List<string>();
                    ml.Add(kp.Key+": Testing begins.");
                    Dictionary<int, Vector.EdgeFamily> uniqueedges = new Dictionary<int, Vector.EdgeFamily>();
                    
                    foreach (SurfaceDefinitions surface in kp.Value)
                    {
                        uniqueedges = Vector.GetEdgeFamilies(surface.SurfaceId, uniqueedges, surface.PlCoords,.0001,.0001);
                    }
                    ml.Add("Gathered edges and their neighboring relationships.");
                    ml.Add("Validating the surfaces' edges alignment with one another (water tightness check.");
                    //check the edge families to see how water tight the edges are
                    //there should always be at least one related Edge
                    //if there is only one, it should match exactly (or within some settable tolerance)
                    //if there is more than one, they each should only intersect at their ends (within some tolerance) and not intersect
                    //and there should be no gaps 
                    //new function added April 11, 2014
                    report = MatchEdges(uniqueedges, ml, report, kp.Key);
                    
                    report.MessageList[kp.Key] = ml;
                }
            }
            catch (Exception e)
            {

            }
            return report;
        }

        public static void FindMatchingEdges(List<gbXMLSpaces.SpaceBoundary>sblist)
        {

            Dictionary<int, DOEgbXMLBasics.EdgeFamily> uniqueedges = new Dictionary<int, DOEgbXMLBasics.EdgeFamily>();
            int distinctedges = 0;
            foreach (gbXMLSpaces.SpaceBoundary sb in sblist)
            {
                int coordcount = sb.sbplane.pl.plcoords.Count;
                for (int i = 0; i < coordcount; i++)
                {
                    //initialize the edge being tested, the test edge
                    DOEgbXMLBasics.EdgeFamily edge = new DOEgbXMLBasics.EdgeFamily();
                    edge.sbdec = sb.surfaceIdRef;
                    edge.relatedEdges = new List<DOEgbXMLBasics.EdgeFamily>();
                    edge.startendpt = new List<Vector.MemorySafe_CartCoord>();
                    if (uniqueedges.Count == 0)
                    {
                        uniqueedges[distinctedges] = edge;
                        //get the first coord in this set, and the coord next to it
                        edge.startendpt.Add(sb.sbplane.pl.plcoords[i]);
                        edge.startendpt.Add(sb.sbplane.pl.plcoords[i+1]);
                        distinctedges++;
                        continue;

                    }
                    //most edges work the same, in terms of the start and end point, except for the last edge (the else case)
                    if (i < coordcount - 1)
                    {
                        edge.startendpt.Add(sb.sbplane.pl.plcoords[i]);
                        edge.startendpt.Add(sb.sbplane.pl.plcoords[i + 1]);
                        //search through existing edges to try and find a perfect match
                        int edgecount = 0; //keeps track of how many guest edges in the dictionary I've searched through
                        foreach(KeyValuePair<int,DOEgbXMLBasics.EdgeFamily> kp in uniqueedges)
                        {
                            
                            Vector.MemorySafe_CartCoord startpt = kp.Value.startendpt[0];
                            //tolerance needed?
                            if (startpt.X == edge.startendpt[0].X && startpt.Y == edge.startendpt[0].Y && startpt.Z == edge.startendpt[0].Z)
                            {
                                //found at least one perfect coordinate match, try to match the second
                                Vector.MemorySafe_CartCoord endpt = kp.Value.startendpt[1];
                                if (endpt.X == edge.startendpt[1].X && endpt.Y == edge.startendpt[1].Y && endpt.Z == edge.startendpt[1].Z)
                                {
                                    //both match, means the match is perfect, so add it to the related surfaces list
                                    kp.Value.relatedEdges.Add(edge);
                                    break;
                                }
                                else
                                {
                                    //the edge may be unique, though it could still have neighboring relationships
                                    //draw vector A
                                    double Ax = endpt.X - edge.startendpt[1].X;
                                    double Ay = endpt.Y - edge.startendpt[1].Y;
                                    double Az = endpt.Z - edge.startendpt[1].Z;
                                    Vector.MemorySafe_CartVect A = new Vector.MemorySafe_CartVect(Ax, Ay, Az);
                                    double Amag = Vector.VectorMagnitude(A);

                                    //take cross product to see if they are even in same plane
                                    double evX = endpt.X - startpt.X;
                                    double evY = endpt.Y - startpt.Y;
                                    double evZ = endpt.Z - startpt.Z;
                                    Vector.MemorySafe_CartVect ev = new Vector.MemorySafe_CartVect(evX, evY, evZ);
                                    double evmag = Vector.VectorMagnitude(ev);
                                    Vector.MemorySafe_CartVect cross = Vector.CrossProduct(A,ev);
                                    double crossmag = Vector.VectorMagnitude(cross);
                                    //tolerance?
                                    if (crossmag == 0)
                                    {
                                        //then we are at least parallel or antiparallel, now see if the point resides on the edge or outside of it
                                        double Bx = startpt.X - edge.startendpt[1].X;
                                        double By = startpt.Y - edge.startendpt[1].Y;
                                        double Bz = startpt.Z - edge.startendpt[1].Z;
                                        Vector.MemorySafe_CartVect B = new Vector.MemorySafe_CartVect(Bx,By,Bz);
                                        double Bmag = Vector.VectorMagnitude(B);
                                        //check to see if the test edge is inside the guest edge
                                        if (Amag < evmag && Bmag < evmag)
                                        {
                                            //this means it lies on the plane at least, so it shares, but it is also still independent because a perfect match wasn't found
                                            kp.Value.relatedEdges.Add(edge);
                                            //accumulate its own relationships
                                            edge.relatedEdges.Add(kp.Value);
                                            edgecount++;
                                            continue;
                                        }

                                        double edgeX = edge.startendpt[1].X - edge.startendpt[0].X;
                                        double edgeY = edge.startendpt[1].Y - edge.startendpt[0].Y;
                                        double edgeZ = edge.startendpt[1].Z - edge.startendpt[0].Z;
                                        Vector.MemorySafe_CartVect edgevec = new Vector.MemorySafe_CartVect(edgeX, edgeY, edgeZ);
                                        double edgemag = Vector.VectorMagnitude(edgevec);

                                        double Cx = startpt.X - edge.startendpt[1].X;
                                        double Cy = startpt.Y - edge.startendpt[1].Y;
                                        double Cz = startpt.Z - edge.startendpt[1].Z;
                                        Vector.MemorySafe_CartVect C = new Vector.MemorySafe_CartVect(Cx, Cy, Cz);
                                        double Cmag = Vector.VectorMagnitude(C);

                                        double Dx = endpt.X - edge.startendpt[1].X;
                                        double Dy = endpt.Y - edge.startendpt[1].Y;
                                        double Dz = endpt.Z - edge.startendpt[1].Z;
                                        Vector.MemorySafe_CartVect D = new Vector.MemorySafe_CartVect(Dx, Dy, Dz);
                                        double Dmag = Vector.VectorMagnitude(D);

                                        if (Dmag < edgemag && Cmag < edgemag)
                                        {
                                            //this means the test edge is longer than the guest edge, but they overlap
                                            kp.Value.relatedEdges.Add(edge);
                                            //the edge is still unique but accumulates a neighbor
                                            edge.relatedEdges.Add(kp.Value);
                                            edgecount++;
                                            continue;
                                        }

                                    }
                                    else
                                    {
                                        //this other point isn't relevant, and the edges don't coincide
                                        edgecount++;
                                        continue;
                                    }

                                }


                            }
                            else if (startpt.X == edge.startendpt[1].X && startpt.Y == edge.startendpt[1].Y && startpt.Z == edge.startendpt[1].Z)
                            {
                                //found at least one perfect coordinate match, try to match the second
                                Vector.MemorySafe_CartCoord endpt = kp.Value.startendpt[1];
                                if (endpt.X == edge.startendpt[0].X && endpt.Y == edge.startendpt[0].Y && endpt.Z == edge.startendpt[0].Z)
                                {
                                    //both match, means the match is perfect, so add it to the related surfaces list
                                    kp.Value.relatedEdges.Add(edge);
                                    break;

                                }
                                else
                                {
                                    //the edge may be unique, though it could still have neighboring relationships
                                    double Ax = endpt.X - edge.startendpt[0].X;
                                    double Ay = endpt.Y - edge.startendpt[0].Y;
                                    double Az = endpt.Z - edge.startendpt[0].Z;
                                    Vector.MemorySafe_CartVect A = new Vector.MemorySafe_CartVect(Ax, Ay, Az);
                                    double Amag = Vector.VectorMagnitude(A);

                                    //take cross product to see if they are even in same plane
                                    double evX = endpt.X - startpt.X;
                                    double evY = endpt.Y - startpt.Y;
                                    double evZ = endpt.Z - startpt.Z;
                                    Vector.MemorySafe_CartVect ev = new Vector.MemorySafe_CartVect(evX, evY, evZ);
                                    double evmag = Vector.VectorMagnitude(ev);
                                    Vector.MemorySafe_CartVect cross = Vector.CrossProduct(A, ev);
                                    double crossmag = Vector.VectorMagnitude(cross);
                                    //tolerance?
                                    if (crossmag == 0)
                                    {
                                        //then we are at least parallel or antiparallel, now see if the point resides on the edge or outside of it
                                        double Bx = startpt.X - edge.startendpt[0].X;
                                        double By = startpt.Y - edge.startendpt[0].Y;
                                        double Bz = startpt.Z - edge.startendpt[0].Z;
                                        Vector.MemorySafe_CartVect B = new Vector.MemorySafe_CartVect(Bx, By, Bz);
                                        double Bmag = Vector.VectorMagnitude(B);
                                        //check to see if the test edge is inside the guest edge
                                        if (Amag < evmag && Bmag < evmag)
                                        {
                                            //this means it lies on the plane at least, so it shares, but it is also still independent because a perfect match wasn't found
                                            kp.Value.relatedEdges.Add(edge);
                                            //accumulate its own relationships
                                            edge.relatedEdges.Add(kp.Value);
                                            edgecount++;
                                            continue;
                                        }

                                        double edgeX = edge.startendpt[1].X - edge.startendpt[0].X;
                                        double edgeY = edge.startendpt[1].Y - edge.startendpt[0].Y;
                                        double edgeZ = edge.startendpt[1].Z - edge.startendpt[0].Z;
                                        Vector.MemorySafe_CartVect edgevec = new Vector.MemorySafe_CartVect(edgeX, edgeY, edgeZ);
                                        double edgemag = Vector.VectorMagnitude(edgevec);

                                        double Cx = startpt.X - edge.startendpt[0].X;
                                        double Cy = startpt.Y - edge.startendpt[0].Y;
                                        double Cz = startpt.Z - edge.startendpt[0].Z;
                                        Vector.MemorySafe_CartVect C = new Vector.MemorySafe_CartVect(Cx, Cy, Cz);
                                        double Cmag = Vector.VectorMagnitude(C);

                                        double Dx = endpt.X - edge.startendpt[0].X;
                                        double Dy = endpt.Y - edge.startendpt[0].Y;
                                        double Dz = endpt.Z - edge.startendpt[0].Z;
                                        Vector.MemorySafe_CartVect D = new Vector.MemorySafe_CartVect(Dx, Dy, Dz);
                                        double Dmag = Vector.VectorMagnitude(D);

                                        if (Dmag < edgemag && Cmag < edgemag)
                                        {
                                            //this means the test edge is longer than the guest edge, but they overlap
                                            kp.Value.relatedEdges.Add(edge);
                                            //the edge is still unique but accumulates a neighbor
                                            edge.relatedEdges.Add(kp.Value);
                                            edgecount++;
                                            continue;
                                        }
                                    }
                                    else
                                    {
                                        //this other point isn't relevant, and the edges don't coincide
                                        edgecount++;
                                        continue;
                                    }
                                }

                            }
                            //neither points perfectly coincide, so we do an exhaustive overlap check.
                            else
                            {
                                Vector.MemorySafe_CartCoord endpt = kp.Value.startendpt[1];
                                //are the two vectors even parallel?  because if they are not, no need to get more complex
                                double evX = endpt.X - startpt.X;
                                double evY = endpt.Y - startpt.Y;
                                double evZ = endpt.Z - startpt.Z;
                                Vector.MemorySafe_CartVect ev = new Vector.MemorySafe_CartVect(evX, evY, evZ);
                                double edgeX = edge.startendpt[1].X - edge.startendpt[0].X;
                                double edgeY = edge.startendpt[1].Y - edge.startendpt[0].Y;
                                double edgeZ = edge.startendpt[1].Z - edge.startendpt[0].Z;
                                Vector.MemorySafe_CartVect edgev = new Vector.MemorySafe_CartVect(edgeX, edgeY, edgeZ);
                                if (Vector.VectorMagnitude(Vector.CrossProduct(ev, edgev)) != 0)
                                {
                                    //they are not even parallel so move on
                                    edgecount++;
                                    continue;
                                }

                                //try to determine if the two edges are parallel
                                //test edge point 1
                                double Ax = endpt.X - edge.startendpt[0].X;
                                double Ay = endpt.Y - edge.startendpt[0].Y;
                                double Az = endpt.Z - edge.startendpt[0].Z;
                                Vector.MemorySafe_CartVect A = new Vector.MemorySafe_CartVect(Ax, Ay, Az);
                                double Amag = Vector.VectorMagnitude(A);

                                //take cross product to see if they are even in same plane
                                evX = endpt.X - startpt.X;
                                evY = endpt.Y - startpt.Y;
                                evZ = endpt.Z - startpt.Z;
                                Vector.MemorySafe_CartVect ev1 = new Vector.MemorySafe_CartVect(evX, evY, evZ);
                                double guestmag = Vector.VectorMagnitude(ev1);
                                Vector.MemorySafe_CartVect cross1 = Vector.CrossProduct(A, ev1);
                                double crossmag = Vector.VectorMagnitude(cross1);
                                //tolerance?
                                if (crossmag == 0)
                                {
                                    //we are at least parallel, now to check for a real intersection
                                    double Bx = startpt.X - edge.startendpt[0].X;
                                    double By = startpt.Y - edge.startendpt[0].Y;
                                    double Bz = startpt.Z - edge.startendpt[0].Z;
                                    Vector.MemorySafe_CartVect B = new Vector.MemorySafe_CartVect(Bx, By, Bz);
                                    double Bmag = Vector.VectorMagnitude(B);
                                    //check to see if the test edge's first point (index 0) is totally inside the guest edge
                                    if (Amag < guestmag && Bmag < guestmag)
                                    {
                                        //the start point of the test edge is inside the guest edge
                                        //test edge point 2 against guest edge point 2
                                        double Cx = endpt.X - edge.startendpt[1].X;
                                        double Cy = endpt.Y - edge.startendpt[1].Y;
                                        double Cz = endpt.Z - edge.startendpt[1].Z;
                                        Vector.MemorySafe_CartVect C = new Vector.MemorySafe_CartVect(Cx, Cy, Cz);
                                        double Cmag = Vector.VectorMagnitude(C);
                                        Vector.MemorySafe_CartVect cross2 = Vector.CrossProduct(C, ev);
                                        crossmag = Vector.VectorMagnitude(cross2);
                                        if (crossmag == 0)
                                        {
                                            //we are at least parallel, in fact we have proven we are totall parallel, now intersect
                                            double Dx = startpt.X - edge.startendpt[1].X;
                                            double Dy = startpt.Y - edge.startendpt[1].Y;
                                            double Dz = startpt.Z - edge.startendpt[1].Z;
                                            Vector.MemorySafe_CartVect D = new Vector.MemorySafe_CartVect(Dx, Dy, Dz);
                                            double Dmag = Vector.VectorMagnitude(D);
                                            if (Cmag < guestmag && Dmag < guestmag)
                                            {
                                                //then it is inside as well, and test vector is engulfed by guest vector
                                                kp.Value.relatedEdges.Add(edge);
                                                //but the edge is still itself unique
                                                edge.relatedEdges.Add(kp.Value);
                                                edgecount++;
                                                continue;
                                            }
                                            else
                                            {
                                                //I am pretty sure that by default, they are still neighbors and this is no difference
                                                //it simply extends beyond one of the ends of the guest vector
                                                kp.Value.relatedEdges.Add(edge);
                                                //but the edge is still itself unique
                                                edge.relatedEdges.Add(kp.Value);
                                                edgecount++;
                                                continue;
                                            }


                                        }
                                        else
                                        {
                                            //we are not parallel, so this is not an adjacency match
                                            edgecount++;
                                            continue;
                                        }
                                    }
                                    else
                                    {
                                        //if test edge start point [index 0] is outside, is one of the guest points inside?
                                        //already computed B
                                        double Cx = startpt.X - edge.startendpt[1].X;
                                        double Cy = startpt.Y - edge.startendpt[1].Y;
                                        double Cz = startpt.Z - edge.startendpt[1].Z;
                                        Vector.MemorySafe_CartVect C = new Vector.MemorySafe_CartVect(Cx, Cy, Cz);
                                        double Cmag = Vector.VectorMagnitude(C);

                                        edgeX = edge.startendpt[1].X - edge.startendpt[0].X;
                                        edgeY = edge.startendpt[1].Y - edge.startendpt[0].Y;
                                        edgeZ = edge.startendpt[1].Z - edge.startendpt[0].Z;
                                        Vector.MemorySafe_CartVect edgevec = new Vector.MemorySafe_CartVect(edgeX, edgeY, edgeZ);
                                        double edgemag = Vector.VectorMagnitude(edgevec);

                                        if (Cmag < edgemag && Bmag < edgemag)
                                        {
                                            //the guest edge's start point is inside the test edge
                                            //guest edge point 2 
                                            double Dx = endpt.X - edge.startendpt[1].X;
                                            double Dy = endpt.Y - edge.startendpt[1].Y;
                                            double Dz = endpt.Z - edge.startendpt[1].Z;
                                            Vector.MemorySafe_CartVect D = new Vector.MemorySafe_CartVect(Dx,Dy,Dz);
                                            double Dmag = Vector.VectorMagnitude(D);
                                            Vector.MemorySafe_CartVect cross3 = Vector.CrossProduct(D, edgevec);
                                            crossmag = Vector.VectorMagnitude(cross3);
                                            if (crossmag == 0)
                                            {
                                                //then we know the two edges are totall parallel and lined up
                                                //determine if the guest edge point 2 is inside the test edge or outside of it
                                                double Ex = startpt.X - edge.startendpt[1].X;
                                                double Ey = startpt.Y - edge.startendpt[1].Y;
                                                double Ez = startpt.Z - edge.startendpt[1].Z;
                                                Vector.MemorySafe_CartVect E = new Vector.MemorySafe_CartVect(Ex, Ey, Ez);
                                                double Emag = Vector.VectorMagnitude(E);
                                                if (Dmag < edgemag && Emag < edgemag)
                                                {
                                                    //it is inside
                                                    kp.Value.relatedEdges.Add(edge);
                                                    //but the edge is still itself unique
                                                    edge.relatedEdges.Add(kp.Value);
                                                    edgecount++;
                                                    continue;
                                                }
                                                else
                                                {
                                                    //it is outside 
                                                    kp.Value.relatedEdges.Add(edge);
                                                    //but the edge is still itself unique
                                                    edge.relatedEdges.Add(kp.Value);
                                                    edgecount++;
                                                    continue;
                                                }
                                            }
                                            else
                                            {
                                                //we are not parallel, so this is not an adjacency match
                                                edgecount++;
                                                continue;
                                            }

                                        }
                                    }



                                }
                                else
                                {
                                    //they are not even parallel, so it is likely best just to shove on
                                    edgecount++;
                                    continue;
                                }


                            }
                        }
                        //this determines if it found a matching edge
                        if (edgecount == uniqueedges.Count)
                        {
                            uniqueedges.Add(distinctedges, edge);
                            distinctedges++;
                        }

                    }
                    //last edge end edge is the zero index   
                    else
                    {
                        edge.startendpt.Add(sb.sbplane.pl.plcoords[i]);
                        edge.startendpt.Add(sb.sbplane.pl.plcoords[0]);
                        int edgecount = 0; //keeps track of how many guest edges in the dictionary I've searched through
                        foreach (KeyValuePair<int, DOEgbXMLBasics.EdgeFamily> kp in uniqueedges)
                        {

                            Vector.MemorySafe_CartCoord startpt = kp.Value.startendpt[0];
                            //tolerance needed?
                            if (startpt.X == edge.startendpt[0].X && startpt.Y == edge.startendpt[0].Y && startpt.Z == edge.startendpt[0].Z)
                            {
                                //found at least one perfect coordinate match, try to match the second
                                Vector.MemorySafe_CartCoord endpt = kp.Value.startendpt[1];
                                if (endpt.X == edge.startendpt[1].X && endpt.Y == edge.startendpt[1].Y && endpt.Z == edge.startendpt[1].Z)
                                {
                                    //both match, means the match is perfect, so add it to the related surfaces list
                                    kp.Value.relatedEdges.Add(edge);
                                    break;
                                }
                                else
                                {
                                    //the edge may be unique, though it could still have neighboring relationships
                                    //draw vector A
                                    double Ax = endpt.X - edge.startendpt[1].X;
                                    double Ay = endpt.Y - edge.startendpt[1].Y;
                                    double Az = endpt.Z - edge.startendpt[1].Z;
                                    Vector.MemorySafe_CartVect A = new Vector.MemorySafe_CartVect(Ax, Ay, Az);
                                    double Amag = Vector.VectorMagnitude(A);

                                    //take cross product to see if they are even in same plane
                                    double evX = endpt.X - startpt.X;
                                    double evY = endpt.Y - startpt.Y;
                                    double evZ = endpt.Z - startpt.Z;
                                    Vector.MemorySafe_CartVect ev = new Vector.MemorySafe_CartVect(evX, evY, evZ);
                                    double evmag = Vector.VectorMagnitude(ev);
                                    Vector.MemorySafe_CartVect cross = Vector.CrossProduct(A, ev);
                                    double crossmag = Vector.VectorMagnitude(cross);
                                    //tolerance?
                                    if (crossmag == 0)
                                    {
                                        //then we are at least parallel or antiparallel, now see if the point resides on the edge or outside of it
                                        double Bx = startpt.X - edge.startendpt[1].X;
                                        double By = startpt.Y - edge.startendpt[1].Y;
                                        double Bz = startpt.Z - edge.startendpt[1].Z;
                                        Vector.MemorySafe_CartVect B = new Vector.MemorySafe_CartVect(Bx, By, Bz);
                                        double Bmag = Vector.VectorMagnitude(B);
                                        //check to see if the test edge is inside the guest edge
                                        if (Amag < evmag && Bmag < evmag)
                                        {
                                            //this means it lies on the plane at least, so it shares, but it is also still independent because a perfect match wasn't found
                                            kp.Value.relatedEdges.Add(edge);
                                            //accumulate its own relationships
                                            edge.relatedEdges.Add(kp.Value);
                                            edgecount++;
                                            continue;
                                        }

                                        double edgeX = edge.startendpt[1].X - edge.startendpt[0].X;
                                        double edgeY = edge.startendpt[1].Y - edge.startendpt[0].Y;
                                        double edgeZ = edge.startendpt[1].Z - edge.startendpt[0].Z;
                                        Vector.MemorySafe_CartVect edgevec = new Vector.MemorySafe_CartVect(edgeX, edgeY, edgeZ);
                                        double edgemag = Vector.VectorMagnitude(edgevec);

                                        double Cx = startpt.X - edge.startendpt[1].X;
                                        double Cy = startpt.Y - edge.startendpt[1].Y;
                                        double Cz = startpt.Z - edge.startendpt[1].Z;
                                        Vector.MemorySafe_CartVect C = new Vector.MemorySafe_CartVect(Cx, Cy, Cz);
                                        double Cmag = Vector.VectorMagnitude(C);

                                        double Dx = endpt.X - edge.startendpt[1].X;
                                        double Dy = endpt.Y - edge.startendpt[1].Y;
                                        double Dz = endpt.Z - edge.startendpt[1].Z;
                                        Vector.MemorySafe_CartVect D = new Vector.MemorySafe_CartVect(Dx, Dy, Dz);
                                        double Dmag = Vector.VectorMagnitude(D);

                                        if (Dmag < edgemag && Cmag < edgemag)
                                        {
                                            //this means the test edge is longer than the guest edge, but they overlap
                                            kp.Value.relatedEdges.Add(edge);
                                            //the edge is still unique but accumulates a neighbor
                                            edge.relatedEdges.Add(kp.Value);
                                            edgecount++;
                                            continue;
                                        }

                                    }
                                    else
                                    {
                                        //this other point isn't relevant, and the edges don't coincide
                                        edgecount++;
                                        continue;
                                    }

                                }


                            }
                            else if (startpt.X == edge.startendpt[1].X && startpt.Y == edge.startendpt[1].Y && startpt.Z == edge.startendpt[1].Z)
                            {
                                //found at least one perfect coordinate match, try to match the second
                                Vector.MemorySafe_CartCoord endpt = kp.Value.startendpt[1];
                                if (endpt.X == edge.startendpt[0].X && endpt.Y == edge.startendpt[0].Y && endpt.Z == edge.startendpt[0].Z)
                                {
                                    //both match, means the match is perfect, so add it to the related surfaces list
                                    kp.Value.relatedEdges.Add(edge);
                                    break;

                                }
                                else
                                {
                                    //the edge may be unique, though it could still have neighboring relationships
                                    double Ax = endpt.X - edge.startendpt[0].X;
                                    double Ay = endpt.Y - edge.startendpt[0].Y;
                                    double Az = endpt.Z - edge.startendpt[0].Z;
                                    Vector.MemorySafe_CartVect A = new Vector.MemorySafe_CartVect(Ax, Ay, Az);
                                    double Amag = Vector.VectorMagnitude(A);

                                    //take cross product to see if they are even in same plane
                                    double evX = endpt.X - startpt.X;
                                    double evY = endpt.Y - startpt.Y;
                                    double evZ = endpt.Z - startpt.Z;
                                    Vector.MemorySafe_CartVect ev = new Vector.MemorySafe_CartVect(evX, evY, evZ);
                                    double evmag = Vector.VectorMagnitude(ev);
                                    Vector.MemorySafe_CartVect cross = Vector.CrossProduct(A, ev);
                                    double crossmag = Vector.VectorMagnitude(cross);
                                    //tolerance?
                                    if (crossmag == 0)
                                    {
                                        //then we are at least parallel or antiparallel, now see if the point resides on the edge or outside of it
                                        double Bx = startpt.X - edge.startendpt[0].X;
                                        double By = startpt.Y - edge.startendpt[0].Y;
                                        double Bz = startpt.Z - edge.startendpt[0].Z;
                                        Vector.MemorySafe_CartVect B = new Vector.MemorySafe_CartVect(Bx, By, Bz);
                                        double Bmag = Vector.VectorMagnitude(B);
                                        //check to see if the test edge is inside the guest edge
                                        if (Amag < evmag && Bmag < evmag)
                                        {
                                            //this means it lies on the plane at least, so it shares, but it is also still independent because a perfect match wasn't found
                                            kp.Value.relatedEdges.Add(edge);
                                            //accumulate its own relationships
                                            edge.relatedEdges.Add(kp.Value);
                                            edgecount++;
                                            continue;
                                        }

                                        double edgeX = edge.startendpt[1].X - edge.startendpt[0].X;
                                        double edgeY = edge.startendpt[1].Y - edge.startendpt[0].Y;
                                        double edgeZ = edge.startendpt[1].Z - edge.startendpt[0].Z;
                                        Vector.MemorySafe_CartVect edgevec = new Vector.MemorySafe_CartVect(edgeX, edgeY, edgeZ);
                                        double edgemag = Vector.VectorMagnitude(edgevec);

                                        double Cx = startpt.X - edge.startendpt[0].X;
                                        double Cy = startpt.Y - edge.startendpt[0].Y;
                                        double Cz = startpt.Z - edge.startendpt[0].Z;
                                        Vector.MemorySafe_CartVect C = new Vector.MemorySafe_CartVect(Cx, Cy, Cz);
                                        double Cmag = Vector.VectorMagnitude(C);

                                        double Dx = endpt.X - edge.startendpt[0].X;
                                        double Dy = endpt.Y - edge.startendpt[0].Y;
                                        double Dz = endpt.Z - edge.startendpt[0].Z;
                                        Vector.MemorySafe_CartVect D = new Vector.MemorySafe_CartVect(Dx, Dy, Dz);
                                        double Dmag = Vector.VectorMagnitude(D);

                                        if (Dmag < edgemag && Cmag < edgemag)
                                        {
                                            //this means the test edge is longer than the guest edge, but they overlap
                                            kp.Value.relatedEdges.Add(edge);
                                            //the edge is still unique but accumulates a neighbor
                                            edge.relatedEdges.Add(kp.Value);
                                            edgecount++;
                                            continue;
                                        }
                                    }
                                    else
                                    {
                                        //this other point isn't relevant, and the edges don't coincide
                                        edgecount++;
                                        continue;
                                    }
                                }

                            }
                            //neither points perfectly coincide, so we do an exhaustive overlap check.
                            else
                            {
                                Vector.MemorySafe_CartCoord endpt = kp.Value.startendpt[1];
                                //are the two vectors even parallel?  because if they are not, no need to get more complex
                                double evX = endpt.X - startpt.X;
                                double evY = endpt.Y - startpt.Y;
                                double evZ = endpt.Z - startpt.Z;
                                Vector.MemorySafe_CartVect ev = new Vector.MemorySafe_CartVect(evX, evY, evZ);
                                double edgeX = edge.startendpt[1].X - edge.startendpt[0].X;
                                double edgeY = edge.startendpt[1].Y - edge.startendpt[0].Y;
                                double edgeZ = edge.startendpt[1].Z - edge.startendpt[0].Z;
                                Vector.MemorySafe_CartVect edgev = new Vector.MemorySafe_CartVect(edgeX, edgeY, edgeZ);
                                if (Vector.VectorMagnitude(Vector.CrossProduct(ev, edgev)) != 0)
                                {
                                    //they are not even parallel so move on
                                    edgecount++;
                                    continue;
                                }
                                //try to determine if the two edges are parallel
                                
                                //test edge point 1
                                double Ax = endpt.X - edge.startendpt[0].X;
                                double Ay = endpt.Y - edge.startendpt[0].Y;
                                double Az = endpt.Z - edge.startendpt[0].Z;
                                Vector.MemorySafe_CartVect A = new Vector.MemorySafe_CartVect(Ax, Ay, Az);
                                double Amag = Vector.VectorMagnitude(A);

                                //take cross product to see if they are even in same plane
                                evX = endpt.X - startpt.X;
                                evY = endpt.Y - startpt.Y;
                                evZ = endpt.Z - startpt.Z;
                                Vector.MemorySafe_CartVect ev1 = new Vector.MemorySafe_CartVect(evX, evY, evZ);
                                double guestmag = Vector.VectorMagnitude(ev);
                                Vector.MemorySafe_CartVect cross1 = Vector.CrossProduct(A, ev);
                                double crossmag = Vector.VectorMagnitude(cross1);
                                //tolerance?
                                if (crossmag == 0)
                                {
                                    //we are at least parallel, now to check for a real intersection
                                    double Bx = startpt.X - edge.startendpt[0].X;
                                    double By = startpt.Y - edge.startendpt[0].Y;
                                    double Bz = startpt.Z - edge.startendpt[0].Z;
                                    Vector.MemorySafe_CartVect B = new Vector.MemorySafe_CartVect(Bx, By, Bz);
                                    double Bmag = Vector.VectorMagnitude(B);
                                    //check to see if the test edge's first point (index 0) is totally inside the guest edge
                                    if (Amag < guestmag && Bmag < guestmag)
                                    {
                                        //the start point of the test edge is inside the guest edge
                                        //test edge point 2 against guest edge point 2
                                        double Cx = endpt.X - edge.startendpt[1].X;
                                        double Cy = endpt.Y - edge.startendpt[1].Y;
                                        double Cz = endpt.Z - edge.startendpt[1].Z;
                                        Vector.MemorySafe_CartVect C = new Vector.MemorySafe_CartVect(Cx, Cy, Cz);
                                        double Cmag = Vector.VectorMagnitude(C);
                                        Vector.MemorySafe_CartVect cross2 = Vector.CrossProduct(C, ev);
                                        crossmag = Vector.VectorMagnitude(cross2);
                                        if (crossmag == 0)
                                        {
                                            //we are at least parallel, in fact we have proven we are totall parallel, now intersect
                                            double Dx = startpt.X - edge.startendpt[1].X;
                                            double Dy = startpt.Y - edge.startendpt[1].Y;
                                            double Dz = startpt.Z - edge.startendpt[1].Z;
                                            Vector.MemorySafe_CartVect D = new Vector.MemorySafe_CartVect(Dx, Dy, Dz);
                                            double Dmag = Vector.VectorMagnitude(D);
                                            if (Cmag < guestmag && Dmag < guestmag)
                                            {
                                                //then it is inside as well, and test vector is engulfed by guest vector
                                                kp.Value.relatedEdges.Add(edge);
                                                //but the edge is still itself unique
                                                edge.relatedEdges.Add(kp.Value);
                                                edgecount++;
                                                continue;
                                            }
                                            else
                                            {
                                                //I am pretty sure that by default, they are still neighbors and this is no difference
                                                //it simply extends beyond one of the ends of the guest vector
                                                kp.Value.relatedEdges.Add(edge);
                                                //but the edge is still itself unique
                                                edge.relatedEdges.Add(kp.Value);
                                                edgecount++;
                                                continue;
                                            }


                                        }
                                        else
                                        {
                                            //we are not parallel, so this is not an adjacency match
                                            edgecount++;
                                            continue;
                                        }
                                    }
                                    else
                                    {
                                        //if test edge start point [index 0] is outside, is one of the guest points inside?
                                        //already computed B
                                        double Cx = startpt.X - edge.startendpt[1].X;
                                        double Cy = startpt.Y - edge.startendpt[1].Y;
                                        double Cz = startpt.Z - edge.startendpt[1].Z;
                                        Vector.MemorySafe_CartVect C = new Vector.MemorySafe_CartVect(Cx, Cy, Cz);
                                        double Cmag = Vector.VectorMagnitude(C);

                                        edgeX = edge.startendpt[1].X - edge.startendpt[0].X;
                                        edgeY = edge.startendpt[1].Y - edge.startendpt[0].Y;
                                        edgeZ = edge.startendpt[1].Z - edge.startendpt[0].Z;
                                        Vector.MemorySafe_CartVect edgevec = new Vector.MemorySafe_CartVect(edgeX, edgeY, edgeZ);
                                        double edgemag = Vector.VectorMagnitude(edgevec);

                                        if (Cmag < edgemag && Bmag < edgemag)
                                        {
                                            //the guest edge's start point is inside the test edge
                                            //guest edge point 2 
                                            double Dx = endpt.X - edge.startendpt[1].X;
                                            double Dy = endpt.Y - edge.startendpt[1].Y;
                                            double Dz = endpt.Z - edge.startendpt[1].Z;
                                            Vector.MemorySafe_CartVect D = new Vector.MemorySafe_CartVect(Dx, Dy, Dz);
                                            double Dmag = Vector.VectorMagnitude(D);
                                            Vector.MemorySafe_CartVect cross3 = Vector.CrossProduct(D, edgevec);
                                            crossmag = Vector.VectorMagnitude(cross3);
                                            if (crossmag == 0)
                                            {
                                                //then we know the two edges are totall parallel and lined up
                                                //determine if the guest edge point 2 is inside the test edge or outside of it
                                                double Ex = startpt.X - edge.startendpt[1].X;
                                                double Ey = startpt.Y - edge.startendpt[1].Y;
                                                double Ez = startpt.Z - edge.startendpt[1].Z;
                                                Vector.MemorySafe_CartVect E = new Vector.MemorySafe_CartVect(Ex, Ey, Ez);
                                                double Emag = Vector.VectorMagnitude(E);
                                                if (Dmag < edgemag && Emag < edgemag)
                                                {
                                                    //it is inside
                                                    kp.Value.relatedEdges.Add(edge);
                                                    //but the edge is still itself unique
                                                    edge.relatedEdges.Add(kp.Value);
                                                    edgecount++;
                                                    continue;
                                                }
                                                else
                                                {
                                                    //it is outside 
                                                    kp.Value.relatedEdges.Add(edge);
                                                    //but the edge is still itself unique
                                                    edge.relatedEdges.Add(kp.Value);
                                                    edgecount++;
                                                    continue;
                                                }
                                            }
                                            else
                                            {
                                                //we are not parallel, so this is not an adjacency match
                                                edgecount++;
                                                continue;
                                            }

                                        }
                                    }



                                }
                                else
                                {
                                    //they are not even parallel, so it is likely best just to shove on
                                    edgecount++;
                                    continue;
                                }


                            }
                        }
                        //this determines if it found a matching edge
                        if (edgecount == uniqueedges.Count)
                        {
                            uniqueedges.Add(distinctedges, edge);
                            distinctedges++;
                        }
                    }

                }
            }
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
