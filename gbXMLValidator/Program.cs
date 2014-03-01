using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Xml;
using System.IO;
using DOEgbXML;
using VectorMath;
using log4net;

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
            string path = "C:\\gbXML\\4_2_0.xml";
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
            //make a simplified representation of the spaces
            List<DOEgbXML.gbXMLSpaces> spaces = DOEgbXML.gbXMLSpaces.getSimpleSpaces(myxml, nsm);
            List<gbXMLSpaces.SpaceBoundary> sblist = gbXMLSpaces.GetSpaceBoundaryList(myxml, nsm);
            
            //check that all polyloops are in a counterclockwise direction
            report = DOEgbXML.gbXMLSpaces.SpaceSurfacesCCTest(spaces, report);
            //process report
            report.Clear();
            //-check for non-planar objects for all Spaces' polyloops
            report = DOEgbXML.gbXMLSpaces.SpaceSurfacesPlanarTest(spaces, report);
            //process report
            report.Clear();

            //4-check for self-intersecting polygons
            report = DOEgbXML.gbXMLSpaces.SpaceSurfacesSelfIntersectionTest(spaces, report);
            //process report
            report.Clear();
            //Vertex Matching------------------------------------------------
            FindMatchingEdges(sblist);

            ////try to parse out the Space Boundary polyloops
            //Dictionary<Vector.CartCoord,Tuple<List<string>,List<bool>>> sbvertices = GetSBVertices(nsm, myxml);
            ////do all vertices have at least one match?  If yes, PASS and move on, if not, then.
            //VertexListToFile(sbvertices,"SpaceBoundaryCoords.txt");
            ////try to parse out the surfaces into my surface objects
            ////check the vertex files
            //report = gbXMLSpaces.findStraySBVertices(@"C:\Temp\gbXML\SpaceBoundaryCoords.txt", report);
            ////process report
            //report.Clear();

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

            
            List<SurfaceDefinitions> surfaces = DOEgbXML.XMLParser.MakeSurfaceList(myxml, nsm);
            //counter clockwise winding test
            report = SurfaceDefinitions.SurfaceCCTest(surfaces, report);
            //process report
            report.Clear();

            //planar surface test
            report = SurfaceDefinitions.TestSurfacePlanarTest(surfaces, report);
            //process report
            report.Clear();

            //self intersecting polygon test
            report = SurfaceDefinitions.SurfaceSelfIntersectionTest(surfaces, report);
            //process the report
            report.Clear();

            //Does the polyloop right hand rule vector form the proper azimuth and tilt? (with and without a CADModelAzimuth)
            report.tolerance = DOEgbXMLBasics.Tolerances.VectorAngleTolerance;
            report = SurfaceDefinitions.SurfaceTiltAndAzCheck(myxml, nsm, report);
            //process report
            report.Clear();

            //Is the Lower Left Corner properly defined?

            //Vertex Matching
            Dictionary<Vector.CartCoord, Tuple<List<string>, List<bool>>> surfvertices = GetSurfVertices(nsm, myxml);
            VertexListToFile(surfvertices, "SurfacesCoords.txt");
            //check the vertex files
            report = gbXMLSpaces.findStraySBVertices(@"C:\Temp\gbXML\SurfacesCoords.txt", report);
            //process report
            report.Clear();

            //Openings Tests-----------------------------------------------------

            //Shading Devices Tests----------------------------------------------


            

        }

        private static Dictionary<Vector.CartCoord, Tuple<List<string>, List<bool>>> GetSBVertices(XmlNamespaceManager nsm, XmlDocument zexml)
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

        public static void FindMatchingEdges(List<gbXMLSpaces.SpaceBoundary>sblist)
        {

            Dictionary<int, DOEgbXMLBasics.EdgeFamily> edges = new Dictionary<int, DOEgbXMLBasics.EdgeFamily>();
            int distinctedges = 0;
            foreach (gbXMLSpaces.SpaceBoundary sb in sblist)
            {
                int coordcount = sb.sbplane.pl.plcoords.Count;
                for (int i = 0; i < coordcount; i++)
                {
                    //test edge
                    DOEgbXMLBasics.EdgeFamily edge = new DOEgbXMLBasics.EdgeFamily();
                    edge.sbdec = sb.surfaceIdRef;
                    edge.relatedEdges = new List<DOEgbXMLBasics.EdgeFamily>();
                    edge.startendpt = new List<Vector.MemorySafe_CartCoord>();
                    if (edges.Count == 0)
                    {
                        edges[distinctedges] = edge;
                        edge.startendpt.Add(sb.sbplane.pl.plcoords[i]);
                        edge.startendpt.Add(sb.sbplane.pl.plcoords[i+1]);
                        edge.sbdec = sb.surfaceIdRef;

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
                        foreach(KeyValuePair<int,DOEgbXMLBasics.EdgeFamily> kp in edges)
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
                        if (edgecount == edges.Count)
                        {
                            edges.Add(distinctedges, edge);
                            distinctedges++;
                        }

                    }
                    //last edge end edge is the zero index   
                    else
                    {
                        edge.startendpt.Add(sb.sbplane.pl.plcoords[i]);
                        edge.startendpt.Add(sb.sbplane.pl.plcoords[0]);
                        int edgecount = 0; //keeps track of how many guest edges in the dictionary I've searched through
                        foreach (KeyValuePair<int, DOEgbXMLBasics.EdgeFamily> kp in edges)
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
                        if (edgecount == edges.Count)
                        {
                            edges.Add(distinctedges, edge);
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
