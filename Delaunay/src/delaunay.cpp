#include "objects.hpp"
#include "sorting.hpp"
#include "delaunay.hpp"
#include "functions.hpp"

using namespace LibraryObjects;
using namespace LibrarySorting;
using namespace LibraryFunctions;

namespace LibraryDelaunay
{
//static constexpr double tolerance = numeric_limits<double>::epsilon();
static constexpr double tolerance = 1e-12;

///DEF Delaunay
//Trova il primo triangolo della triangolazione
void Delaunay::firstTriangle(vector<Point>& sortedX, vector<Point>& sortedY, vector<Triangle>& triangles, vector<Segment>& convexHull, vector<Point>& hullPoints, vector<int>& hullTrianglesIndices) {

    //Creo i 4 triangoli che si possono ottenere combinando i 4 estremi dei due vettori
    //attenzione:               //se due estremi coincidono ho solo 1 triangolo
                                //se due coppie di estremi coincidono non ho abbastanza punti

    Point p1 = sortedX[0];
    Point p2 = sortedX[sortedX.size()-1];
    Point p3 = sortedY[0];
    Point p4 = sortedY[sortedY.size()-1];
    vector<Point> vectorPoints;
    vectorPoints = {p1, p2, p3, p4};

    vector<Point> distinctPoints;
    //Verifica se dato un gruppo di 4 punti, questi sono diversi tra loro
    if (!areSamePoint(p1,p2) && !areSamePoint(p1,p3) && !areSamePoint(p1,p4)) {
        distinctPoints.push_back(p1);
    }
    if (!areSamePoint(p2,p3) && !areSamePoint(p2,p4)) {
        distinctPoints.push_back(p2);
    }
    if (!areSamePoint(p3,p4)) {
        distinctPoints.push_back(p3);
        distinctPoints.push_back(p4);
    }
    else {
        distinctPoints.push_back(p3);
    }

    //creo un triangolo vuoto che diventerà il triangolo migliore (cioè quello con l'area più grande)
    Triangle bestTriangle;

    //se ho 4 punti distinti:
    if (distinctPoints.size() == 4)
    { vector<Triangle> tryTriangles;

        //costruisco i 4 triangoli candidati con i miei 4 punti:
        Triangle t1 = Triangle(p1, p2, p3);
        tryTriangles.push_back(Triangle(p1, p2, p3));
        tryTriangles.push_back(Triangle(p1, p2, p4));
        tryTriangles.push_back(Triangle(p2, p3, p4));
        tryTriangles.push_back(Triangle(p1, p3, p4));

        //calcolo le aree dei triangoli e mi tengo la piu grande
        double bestArea = 0;
        double areaTriangle;

    for (unsigned int i=0; i<4; i++) {
        areaTriangle = tryTriangles[i].Area(tryTriangles[i].p1, tryTriangles[i].p2, tryTriangles[i].p3);
        if (areaTriangle>bestArea){
            bestArea = areaTriangle;
            bestTriangle = tryTriangles[i];
            }
        }
    }

    //se ho 3 punti distinti, basta costruire il triangolo grande con quei 3 punti
    else if (distinctPoints.size() == 3) {
        bestTriangle = Triangle(distinctPoints[0], distinctPoints[1], distinctPoints[2]);
    }

    //se ho 2 punti distinti:
    else if (distinctPoints.size() == 2) {
        Triangle triangle3;
        Triangle triangle4;
        Point newp3 = sortedX[1];  //prendo il secondo punto più piccolo secondo x
        Point newp4 = sortedX[sortedX.size() - 2];  //prendo il penultimo punto più piccolo secondo x
        double areap3 = triangle3.Area(distinctPoints[0], distinctPoints[1], newp3);  //area con il secondo punto
        double areap4 = triangle4.Area(distinctPoints[0], distinctPoints[1], newp4);  //area con il penultimo punto
        if (areap3<areap4)
            bestTriangle = Triangle(distinctPoints[0], distinctPoints[1], newp4);
        else
            bestTriangle = Triangle(distinctPoints[0], distinctPoints[1], newp3);
    }

    //tolgo da sortedX i 3 vertici del triangolo grande appena creato
    sortedX.erase(sortedX.begin() + binarySearch(sortedX, bestTriangle.p1));
    sortedX.erase(sortedX.begin() + binarySearch(sortedX, bestTriangle.p2));
    sortedX.erase(sortedX.begin() + binarySearch(sortedX, bestTriangle.p3));

    //aggiungo questo triangolo ordinato in senso antiorario al vettore dei triangoli
    antiClockWiseOrder(bestTriangle);
    triangles.push_back(bestTriangle);

    //agigungo i punti al guscio
    hullPoints.push_back(bestTriangle.p1);
    hullPoints.push_back(bestTriangle.p2);
    hullPoints.push_back(bestTriangle.p3);

    //aggiungo segmenti di bordo
    convexHull.push_back(Segment(bestTriangle.p1, bestTriangle.p2));
    convexHull[0].coefficients = convexHull[0].lineCoefficients(bestTriangle.p1, bestTriangle.p2);

    convexHull.push_back(Segment(bestTriangle.p2, bestTriangle.p3));
    convexHull[1].coefficients = convexHull[1].lineCoefficients(bestTriangle.p2, bestTriangle.p3);

    convexHull.push_back(Segment(bestTriangle.p3, bestTriangle.p1));
    convexHull[2].coefficients = convexHull[2].lineCoefficients(bestTriangle.p3, bestTriangle.p1);

    //neighborhood settato di default
    hullTrianglesIndices.push_back(0);

}

//Verifica se point è interno/esterno al circumcerchio di ogni traingolo.
//Si ferma o quando trova un circumcerchio a cui è interno o quando li ha verificati tutti e sta fuori
int Delaunay::isPointInsideCircumcircle(Point& point, vector<Triangle>& triangles, int& triangleIndex) {

    int triangleSize = triangles.size();
    for (int i = triangleIndex; i<triangleSize; i++) {
        //seleziona triangolo
        Triangle& triangle = triangles[i];

        double ax = triangle.p1.x - point.x; //Ax-Dx
        double ay = triangle.p1.y - point.y; //Ax-Dy
        double bx = triangle.p2.x - point.x; //Bx-Dx
        double by = triangle.p2.y - point.y; //Bx-Dy
        double cx = triangle.p3.x - point.x; //Cx-Dx
        double cy = triangle.p3.y - point.y; //Cx-Dy

        double d = ax * (by * (pow(cx,2) + pow(cy,2)) - cy * (pow(bx,2) + pow(by,2))) -
                   bx * (ay * (pow(cx,2) + pow(cy,2)) - cy * (pow(ax,2) + pow(ay,2))) +
                   cx * (ay * (pow(bx,2) + pow(by,2)) - by * (pow(ax,2) + pow(ay,2)));

        if (d > 0)   //d>0 dice che il punto sta dentro circumucerchio del triangolo
            return i;  //restituisce indice di triangolo appena testato
    }

    return triangleIndex = -1;
}

//Verifica se un punto è interno/esterno al triangolo rilevato dal circumcerchio di isPointInsideCircumcircle
 int Delaunay::isPointInsideTriangle(Point& point, vector<Triangle>& triangles, int& triangleIndex) {

    //prendo il triangolo dal vettore triangles
    Triangle triangle = triangles[triangleIndex];

    //determinante che serve per determinare se il punto è interno al triangolo
    double d1 = (point.x - triangle.p1.x) * (triangle.p3.y - triangle.p1.y) - (point.y - triangle.p1.y) * (triangle.p3.x - triangle.p1.x);
    double d2 = (triangle.p1.x - point.x) * (triangle.p2.y - triangle.p1.y) - (triangle.p1.y - point.y) * (triangle.p2.x - triangle.p1.x);
    double d3 = (triangle.p2.x - triangle.p1.x) * (triangle.p3.y - triangle.p1.y) - (triangle.p2.y - triangle.p1.y) * (triangle.p3.x - triangle.p1.x);
    double a = d1/d3;
    double b = d2/d3;

    //if (a>0 && b>0 && (a+b)<1)    //se tutti 3 determinanti sono positivi, si trova all’interno del triangolo
    if (a > 0 && b > 0 && (a + b) < 1)
        return triangleIndex;       //restituisce triangleIndex (che davamo come input alla funzione) se il punto è interno al triangolo
    else
        return -1;        //altrimenti -1 se è esterno alla triangolazione
 }

 //unisce le funzioni isPointInsideCircumcircle e isPointInsideTriangle
 int Delaunay::findTriangleContainingPoint(Point& point, vector<Triangle>& triangles, int& triangleIndex) {

    if (isPointInsideCircumcircle(point, triangles, triangleIndex)!=-1) {  //se il punto è interno al circumcerchio
        if (isPointInsideTriangle(point, triangles, triangleIndex)!=-1)  //se il punto è interno al triangolo
            return triangleIndex;
        else {  //se il punto è esterno al triangolo ripeto
            int triangleSizeMinus1 = triangles.size()-1;
            if (triangleIndex<triangleSizeMinus1) {
                triangleIndex = triangleIndex + 1;
                findTriangleContainingPoint(point, triangles, triangleIndex);
                return triangleIndex;
            }
            else
                return triangleIndex = -1;
        }
        return triangleIndex;
    }
    else  //se il punto è esterno al circumcerchio
        return triangleIndex = -1;
 }

//Tratta caso punto interno triangolo
//unisce punto interno con 3 vertici x creare sottotriangolazione
 void Delaunay::inTriangle(Point& point, vector<Triangle>& triangles, int& triangleIndex, vector<int>& hullTrianglesIndices){

    //prendo il triangolo dal vettore triangles
    Triangle triangle = triangles[triangleIndex];

    //salvo i punti del triangolo nelle variabili v1, v2, v3
    Point v1 = triangle.p1;
    Point v2 = triangle.p2;
    Point v3 = triangle.p3;

    //creo 3 nuovi triangoli vuoti
    Triangle t1, t2, t3;
    int triangleIndex2 = triangles.size();
    int triangleIndex3 = triangles.size() + 1;

    t1.p1 = v1;    //primo vertice del triangolo grande
    t1.p2 = v2;    //secondo vertice del triangolo grande
    t1.p3 = point;     //il terzo vertice è il punto interno al triangolo grande
    t1.neighbors[0] = triangle.neighbors[0];   //primo lato + l'indice del triangolo vicino del triangolo grande che ha in comune il lato con t1
    t1.neighbors[1] = triangleIndex2;       //secondo lato + indice di t2, che subito dopo creerò
    t1.neighbors[2] = triangleIndex3;       //terzo lato + indice di t3, che subito dopo creerò

    //le adiacenze del "vecchio" triangolo adiacente a t1 rimangono invariate

    t2.p1 = v2;    //secondo vertice del triangolo grande
    t2.p2 = v3;    //terzo vertice del triangolo grande
    t2.p3 = point;     //il terzo vertice è il punto interno al triangolo grande
    t2.neighbors[0] = triangle.neighbors[1];   //primo lato + l'indice del triangolo vicino del triangolo grande che ha in comune il lato con t2
    t2.neighbors[1] = triangleIndex3;       //secondo lato + indice di t3
    t2.neighbors[2] = triangleIndex;           //terzo lato + indice di t1

    //aggiorniamo le adiacenze del triangolo vicino a t2
    for (int i=0; i<3; i++) {
        if (t2.neighbors[0] != -1) {
            if (triangles[t2.neighbors[0]].neighbors[i] == triangleIndex) {
                triangles[t2.neighbors[0]].neighbors[i] = triangleIndex2;
                break;
            }
        }
    }

    t3.p1 = v3;
    t3.p2 = v1;
    t3.p3 = point;
    t3.neighbors[0] = triangle.neighbors[2];  //primo lato + indice del triangolo vicino del triangolo grande che ha in comune il lato con t3
    t3.neighbors[1] = triangleIndex;          //secondo lato + indice di t1
    t3.neighbors[2] = triangleIndex2;      //terzo lato + indice di t2

    //aggiorniamo le adiacenze del triangolo vicino a t3
    for (int i=0; i<3; i++) {
        if (t3.neighbors[0] != -1) {
            if (triangles[t3.neighbors[0]].neighbors[i] == triangleIndex) {
                triangles[t3.neighbors[0]].neighbors[i] = triangleIndex3;
                break;
            }
        }
    }

    //inserisco la sottotriangolazione t1 t2 e t3 dentro il vettore triangles
    triangles[triangleIndex] = t1;  //inserisco t1 al posto del triangolo grande
    triangles.push_back(t2);  //inserisco t2 alla fine
    triangles.push_back(t3);  //inserisco t3 alla fine

    /*se uno dei 3 nuovi triangoli non ha vicini, significa che sta sul guscio: rimuoviamo dunque triangleIndex dal vettore contenente
    gli indici dei triangoli sul guscio*/
    if (t1.neighbors[0]==-1 || t2.neighbors[0]==-1 || t3.neighbors[0]==-1) {
        hullTrianglesIndices.erase(remove(hullTrianglesIndices.begin(), hullTrianglesIndices.end(), triangleIndex), hullTrianglesIndices.end());
        if (t1.neighbors[0]==-1)  //se t1 non ha vicini, inserisco il suo indice (che adesso è triangleIndex) dentro hullTrianglesIndices
            hullTrianglesIndices.push_back(triangleIndex);
        if (t2.neighbors[0]==-1)  //se t2 non ha vicini, inserisco il suo indice (che corrisponde al penultimo elemento di triangles) dentro hullTrianglesIndices
            hullTrianglesIndices.push_back(triangleIndex2);
        if (t3.neighbors[0]==-1)  //se t3 non ha vicini, inserisco il suo indice (che corrisponde all'ultimo elemento di triangles) dentro hullTrianglesIndices
            hullTrianglesIndices.push_back(triangleIndex3);
    }

    verifyDelaunayCondition(t1, triangleIndex, triangles, hullTrianglesIndices);
    //triangles[triangleIndex] = t1;
    verifyDelaunayCondition(t2, triangleIndex2, triangles, hullTrianglesIndices);
    //triangles[triangleIndex2] = t2;
    verifyDelaunayCondition(t3, triangleIndex3, triangles, hullTrianglesIndices);
    //triangles[triangleIndex3] = t3;
 }

//Caso punto esterno triangolazione
//unisce punto esterno con tutti i vertici su guscio e rimuove segmenti che hanno intersezione
//aggiorna neighbors e guscio
 void Delaunay::outTriangle(Point& outsidePoint, vector<Triangle>& triangles, vector<Segment>& convexHull, vector<Point>& hullPoints, vector<int>& hullTrianglesIndices){


     //creo un vettore che conterrà i nuovi segmenti
     vector<Segment> newSegments;
     Segment temporary;  //creazione di segmento temporaneo dove salvo i coefficienti
     for (unsigned int i=0; i<hullPoints.size(); i++){  //cicliamo su tutti i punti del guscio
          //creiamo una retta (vettore con m e q) dato il punto esterno e l'i-esimo punto del guscio
          Vector2d possibleNewLine = temporary.lineCoefficients(outsidePoint, hullPoints[i]);
          /*iniziamo a contare il numero di intersezioni tra l'i-esima retta (corrispontende all'i-esimo punto)
          e le rette che corrispondono ai segmenti del guscio*/
          int numberOfIntersections = 0;
          for (unsigned int j=0; j<hullPoints.size(); j++){
              //troviamo il punto di intersezione tra l'i-esima retta e la j-esima retta (corrispondente al j-esimo segmento del guscio):
              Point intersection = findIntersection(possibleNewLine, convexHull[j].coefficients);
              //se almeno una delle due coordinate esiste:
              if (!isnan(intersection.x) || !isnan(intersection.y)){
                 //calcolo la lunghezza del j-esimo segmento del guscio
                 double lenghtSegment = getDistance(convexHull[j].p1, convexHull[j].p2);
                 /*se entrambe le distanze tra il punto intersezione e i due punti del j-esimo segmento
                 sono minori della lunghezza del j-esimo segmento, ho intersezione e devo aggiungerla a numberOfIntersections:*/

                 double distance1 = getDistance(intersection, convexHull[j].p1);
                 double distance2 = getDistance(intersection, convexHull[j].p2);

                 //se m e q delle due rette confrontate per intersezione sono uguali, allora le rette sono coincidenti. l'intersezione va ritenuta valida
                 if (areEqual(possibleNewLine, convexHull[j].coefficients))
                     numberOfIntersections++;
                 else if (distance1<=(lenghtSegment+tolerance) && distance2<=(lenghtSegment+tolerance)) {
                     double lengthSide = getDistance(outsidePoint, hullPoints[i]);
                     double distance = getDistance(outsidePoint, intersection);
                     if (lengthSide >= (distance-tolerance))
                         numberOfIntersections++;
                     /*se il numero di intersezioni supera 2, non posso costruire il triangolo usando
                     il punto esterno e i due punti del j-esimo segmento: uscirò dunque dal ciclo con j e passo al prossimo i*/
                 }
                 if (numberOfIntersections>2)
                     break;
                     /*se ho 2 intersezioni, costruisco il nuovo segmento con il punto esterno
                     e l'i-esimo punto del guscio; lo aggiungo a newSegments*/
              }
          }
          if (numberOfIntersections==2)
             newSegments.push_back(Segment(outsidePoint, hullPoints[i], possibleNewLine));

      }

    //aggiungo il nuovo punto a hullPoints
    hullPoints.push_back(outsidePoint);

    int numberOfNewTriangles = 0;
    for (unsigned int i=0; i<newSegments.size(); i++){
       for (unsigned int j=0; j<newSegments.size(); j++){

            if (i<j) {
                //creo un segmento con due punti della convexHull
                Segment verifySegment = Segment(newSegments[i].p2, newSegments[j].p2);
                //Vector2d verifyCoefficients = verifySegment.lineCoefficients(newSegments[i].p2, newSegments[j].p2);
                //verifySegment.coefficients = verifyCoefficients;

                bool existed = false;
                for (unsigned int q=0; q<convexHull.size(); q++) {
                    if (areSameSegment(verifySegment, convexHull[q])) {
                        existed = true;
                        break;
                    }
                 }

                //verifico se segmento che chiude triangolo costruito con i due segmenti che sto accoppiando fa parte del guscio, qui == quindi SI
                if (existed) {
                    Triangle newTriangle = Triangle(outsidePoint, newSegments[i].p2, newSegments[j].p2);
                    antiClockWiseOrder(newTriangle);
                    int newTriangleIndex = triangles.size();
                    triangles.push_back(newTriangle);
                    numberOfNewTriangles = numberOfNewTriangles+1;

                   //aggiorno le adiacenze
                   unsigned int temporary = 0;
                   Triangle hullTriangle;
                   int hullTriangleIndicesSize = hullTrianglesIndices.size();
                   for (int k=0; k<hullTriangleIndicesSize; k++){                                                     //per ogni triangolo del guscio controllo se il segmento di adiacenza con il nuovo triangolo aggiunto gli appartiene (newTriangle.p1 = outsidePoint)
                        hullTriangle = triangles[hullTrianglesIndices[k]];
                        if (areSamePoint(newTriangle.p2, hullTriangle.p1) &&                                        //confronto questi vertici di hullTriangle perchè so che sono ordinati in senso antiorario
                            areSamePoint(newTriangle.p3, hullTriangle.p3))  {
                            hullTriangle.neighbors[2] = newTriangleIndex;                                    //aggiorno l'indice del vicino mettendogli l'indice del triangolo appena aggiunto
                            temporary = k;
                            break;
                        }

                        else if (areSamePoint(newTriangle.p2, hullTriangle.p2) &&                                        //confronto questi vertici di hullTriangle perchè so che sono ordinati in senso antiorario
                            areSamePoint(newTriangle.p3, hullTriangle.p1))  {
                            hullTriangle.neighbors[0] = newTriangleIndex;                                    //aggiorno l'indice del vicino mettendogli l'indice del triangolo appena aggiunto
                            temporary = k;
                            break;
                        }
                        else if (areSamePoint(newTriangle.p2, hullTriangle.p3) &&                                        //confronto questi vertici di hullTriangle perchè so che sono ordinati in senso antiorario
                            areSamePoint(newTriangle.p3, hullTriangle.p2))  {
                            hullTriangle.neighbors[1] = newTriangleIndex;                                    //aggiorno l'indice del vicino mettendogli l'indice del triangolo appena aggiunto
                            temporary = k;
                            break;
                        }
                        //se non era indice k giusto, esco e continuo con il successivo
                    }
                    triangles[hullTrianglesIndices[temporary]].neighbors = hullTriangle.neighbors;
                    newTriangle.neighbors[1] = hullTrianglesIndices[temporary];                                       //ora aggiorno il vicino del nuovo triangolo dal lato adiacente con un vecchio
                    triangles[newTriangleIndex].neighbors = newTriangle.neighbors;

                    //controlliamo se il triangolo adiacente a quello creato ha lati sul bordo
                    bool onSide = false;
                    for (int l = 0; l<3; l++) {
                        if (triangles[hullTrianglesIndices[temporary]].neighbors[l] == -1) {
                            onSide = true;
                            break;
                        }
                    }
                    //se non ha più lati sul bordo, eliminiamo l'indice da hullTrianglesIndices
                    if (onSide == false)
                        hullTrianglesIndices.erase(remove(hullTrianglesIndices.begin(), hullTrianglesIndices.end(), hullTrianglesIndices[temporary]), hullTrianglesIndices.end());      //rimuovo l'indice del triangolo adiacente che prima era di bordo e ora non è più

                    //elimino il lato che prima era nella convexHull e ora non lo è più
                    Segment oldHullSegment = Segment(newTriangle.p2, newTriangle.p3);
                    int convexHullSize = convexHull.size();
                    for (int k=0; k<convexHullSize; k++) {
                        if (areSameSegment(oldHullSegment, convexHull[k])) {
                            convexHull.erase(next(convexHull.begin(), k));
                            break;
                        }
                    }
                }
            }
       }
    }

    //aggiorno i vicini dei nuovi triangoli inseriti solo dai lati di contatto tra i nuovi triangoli
    if (numberOfNewTriangles == 1) {
        int trianglesSizeMinus1 = triangles.size()-1;
        Point p1 = triangles[trianglesSizeMinus1].p1;
        Point p2 = triangles[trianglesSizeMinus1].p2;
        Point p3 = triangles[trianglesSizeMinus1].p3;
        Segment s1 = Segment(p1, p2);
        Segment s2 = Segment(p3, p1);
        s1.coefficients = s1.lineCoefficients(p1, p2);
        s2.coefficients = s2.lineCoefficients(p3, p1);
        convexHull.push_back(s1);
        convexHull.push_back(s2);
    }

    else {
        for (int i=0; i<numberOfNewTriangles; i++) {
            for (int j=0; j<numberOfNewTriangles; j++) {
                if (i!=j) {                                                                                     //così non confronto con se stesso lo stesso triangolo!
                    Triangle triangleToCompare1 = triangles[triangles.size()-numberOfNewTriangles + i];         //salvo temporaneamente per facilitare la lettura del codice
                    Triangle triangleToCompare2 = triangles[triangles.size()-numberOfNewTriangles + j];
                    if (areSamePoint(triangleToCompare1.p2, triangleToCompare2.p3)) {
                        triangleToCompare1.neighbors[0] = triangles.size()-numberOfNewTriangles + j;            //aggiorno in triangles l'indice del triangolo corrispondente a triangleToCompare1 e gli assegno l'indice di quello di fianco a destra (corrispondente a triangleToCompare2)
                        triangleToCompare2.neighbors[2] = triangles.size()-numberOfNewTriangles + i;

                        int hullPointsSize = hullPoints.size();
                        for (int k=0; k<hullPointsSize; k++) {
                            //elimino il punto in comune fra i due triangoli creati perchè non è più di bordo
                            if (areSamePoint(triangleToCompare1.p2, hullPoints[k])) {
                                hullPoints.erase(next(hullPoints.begin(), k));
                                break;
                            }
                        }

                        if (triangleToCompare1.neighbors[2]==-1) {                                                //aggiungo segmento di bordo a convexHull
                            Segment convexHullSegment = Segment(triangleToCompare1.p3, outsidePoint);
                            convexHullSegment.coefficients = convexHullSegment.lineCoefficients(triangleToCompare1.p3 ,outsidePoint);
                            convexHull.push_back(convexHullSegment);
                        }
                        if (triangleToCompare1.neighbors[0]==-1) {
                            Segment convexHullSegment = Segment(outsidePoint, triangleToCompare1.p2);
                            convexHullSegment.coefficients = convexHullSegment.lineCoefficients(outsidePoint, triangleToCompare1.p2);
                            convexHull.push_back(convexHullSegment);
                        }

                        if (triangleToCompare2.neighbors[2]==-1) {                                                //aggiungo segmento di bordo a convexHull
                            Segment convexHullSegment = Segment(triangleToCompare2.p3, outsidePoint);
                            convexHullSegment.coefficients = convexHullSegment.lineCoefficients(triangleToCompare2.p3 ,outsidePoint);
                            convexHull.push_back(convexHullSegment);
                        }
                        if (triangleToCompare2.neighbors[0]==-1) {
                            Segment convexHullSegment = Segment(outsidePoint, triangleToCompare2.p2);
                            convexHullSegment.coefficients = convexHullSegment.lineCoefficients(outsidePoint, triangleToCompare2.p2);
                            convexHull.push_back(convexHullSegment);
                        }

                        triangles[triangles.size()-numberOfNewTriangles + i].neighbors[0] = triangleToCompare1.neighbors[0];
                        triangles[triangles.size()-numberOfNewTriangles + j].neighbors[2] = triangleToCompare2.neighbors[2];
                    }
                }
            }
        }
    }

    int firstNewTriangleIndex = triangles.size()-numberOfNewTriangles;
    for (int i=0; i<numberOfNewTriangles; i++) {
        int triangleIndex = firstNewTriangleIndex+i;
        for (int j=0; j<3; j++) {
            if (triangles[triangleIndex].neighbors[j] == -1) {
                hullTrianglesIndices.push_back(triangleIndex);
                break;
            }
        }
    }


    for (int i=0; i<numberOfNewTriangles; i++) {
        int triangleIndex = firstNewTriangleIndex+i;
        verifyDelaunayCondition(triangles[triangleIndex], triangleIndex, triangles, hullTrianglesIndices);
    }
}

 //Verifico la condizione di Delunay
 void Delaunay::verifyDelaunayCondition(Triangle& triangle, int& triangleIndex, vector<Triangle>& triangles, vector<int>& hullTrianglesIndices){

    vector<int> flippedTrianglesIndices;

    double angle1;
    double angle2;
    Triangle nearTriangle;
    //i sarà l'indice che identifica il lato del triangolo 1 su cui ho verificato avere adiacenza con triangolo 2
    //j sarà l'indice che identifica il lato del triangolo 2 su cui ho verificato avere adiacenza con triangolo 1
    for (int i=0; i<3; i++) {                    //per ogni lato
        if (triangle.neighbors[i]!=-1) {      //triangolo testato ha un triangolo adiacente su quel lato
            nearTriangle = triangles[triangle.neighbors[i]];             //il triangolo vicino è quello indicato dall'indice contenuto in neighbors
            for (int j=0; j<3; j++) {
                if (nearTriangle.neighbors[j] == triangleIndex) {         //capisco su quale lato del triangolo vicino i due triangoli sono adiacenti
                    if (i==0) {
                        Segment segment1 = Segment(triangle.p1, triangle.p3);
                        Segment segment2 = Segment(triangle.p2, triangle.p3);
                        angle1 = getAngle(segment1,segment2);
                    }
                        //angle1 = getAngle(Segment(triangle.p1, triangle.p3), Segment(triangle.p2, triangle.p3));                   //calcolo l'angolo opposto al lato adiacente nel triangolo di partenza
                    else if (i==1)
                        angle1 = getAngle(Segment(triangle.p2, triangle.p1), Segment(triangle.p3, triangle.p1));
                    else if (i==2)
                        angle1 = getAngle(Segment(triangle.p3, triangle.p2), Segment(triangle.p1, triangle.p2));

                    if (j==0)
                        angle2 = getAngle(Segment(nearTriangle.p1, nearTriangle.p3), Segment(nearTriangle.p2, nearTriangle.p3));   //calcolo l'angolo opposto al lato adiacente nel triangolo vicino
                    else if (j==1)
                        angle2 = getAngle(Segment(nearTriangle.p2, nearTriangle.p1), Segment(nearTriangle.p3, nearTriangle.p1));
                    else if (j==2)
                        angle2 = getAngle(Segment(nearTriangle.p3, nearTriangle.p2), Segment(nearTriangle.p1, nearTriangle.p2));

                    //se somma è maggiore di 180° faccio flip
                    if (angle1 + angle2 > 180.0) {
                        int triangle2Index = triangle.neighbors[i];
                        flippedTrianglesIndices.push_back(triangle2Index);
                        flipTriangles(triangleIndex, triangle2Index, i, j, hullTrianglesIndices, triangles);
                        triangle = triangles[triangleIndex];
                        i = 0;
                        break;
                    }

                }
            }
        }
    }
    int numberOfFlips = flippedTrianglesIndices.size();
    for (int k=0; k<numberOfFlips; k++)
        verifyDelaunayCondition(triangles[flippedTrianglesIndices[k]], flippedTrianglesIndices[k], triangles, hullTrianglesIndices);

 }


//Flippo sul triangolo
//l'ordine con cui passo gli indici dei triangoli da flippare è ininfluente
 void Delaunay::flipTriangles(int& triangle1Index, int& triangle2Index, int& adjacentSide1, int& adjacentSide2,
                              vector<int>& hullTrianglesIndices,
                              vector<Triangle>& triangles){

    Triangle oldTriangle1 = triangles[triangle1Index];
    Triangle oldTriangle2 = triangles[triangle2Index];
    //salvo quanti lati sul bordo avevano triangle1 e triangle2 prima di essere flippati. Mi servirà per dopo
    int nEdgeOLDTriangle1 = 0;
    int nEdgeOLDTriangle2 = 0;
    for (int k=0; k<3; k++) {
        if (oldTriangle1.neighbors[k]==-1)
            nEdgeOLDTriangle1++;
        if (oldTriangle2.neighbors[k]==-1)
            nEdgeOLDTriangle2++;
    }

    //in base al lato che ha generato adiacenza salvo i valori dei punti
    Point v1Common;
    Point v2Common;
    Point v3FromTriangle1;
    Point v3FromTriangle2;

    if (adjacentSide1==0) {
        v1Common = oldTriangle1.p1;
        v2Common = oldTriangle1.p2;
        v3FromTriangle1 = oldTriangle1.p3;
    }
    else if (adjacentSide1==1) {
        v1Common = oldTriangle1.p2;
        v2Common = oldTriangle1.p3;
        v3FromTriangle1 = oldTriangle1.p1;
    }
    else if (adjacentSide1==2) {
        v1Common = oldTriangle1.p3;
        v2Common = oldTriangle1.p1;
        v3FromTriangle1 = oldTriangle1.p2;
    }

    if (adjacentSide2==0)
        v3FromTriangle2 = oldTriangle2.p3;
    else if (adjacentSide2==1)
        v3FromTriangle2 = oldTriangle2.p1;
    else if (adjacentSide2==2)
        v3FromTriangle2 = oldTriangle2.p2;

    //devo salvare prima i vicini, poi sostituisco triangolo e poi ripristino vicini
    vector<int> neighborsTriangle1 = oldTriangle1.neighbors;
    vector<int> neighborsTriangle2 = oldTriangle2.neighbors;
    Triangle triangle1 = Triangle(v3FromTriangle1, v1Common, v3FromTriangle2);       //dovrebbero già essere ordinati in senso antiorario
    Triangle triangle2 = Triangle(v3FromTriangle2, v2Common, v3FromTriangle1);

    if (adjacentSide1==0) {
        triangle1.neighbors[0] = neighborsTriangle1[2];
        triangle2.neighbors[1] = neighborsTriangle1[1];

        //aggiorniamo le adiacenze dei triangoli vicini a triangle1
        for (int i=0; i<3; i++) {
            if (neighborsTriangle1[1] != -1) {
                if (triangles[neighborsTriangle1[1]].neighbors[i] == triangle1Index) {
                    triangles[neighborsTriangle1[1]].neighbors[i] = triangle2Index;
                    break;
                }
            }
        }
    }

    else if (adjacentSide1==1) {
        triangle1.neighbors[0] = neighborsTriangle1[0];
        triangle2.neighbors[1] = neighborsTriangle1[2];

        //aggiorniamo le adiacenze dei triangoli vicini a triangle1
        for (int i=0; i<3; i++) {
            if (neighborsTriangle1[2] != -1) {
                if (triangles[neighborsTriangle1[2]].neighbors[i] == triangle1Index) {
                    triangles[neighborsTriangle1[2]].neighbors[i] = triangle2Index;
                    break;
                }
            }
        }
    }

    else if (adjacentSide1==2) {
        triangle1.neighbors[0] = neighborsTriangle1[1];
        triangle2.neighbors[1] = neighborsTriangle1[0];

        //aggiorniamo le adiacenze dei triangoli vicini a triangle1
        for (int i=0; i<3; i++) {
            if (neighborsTriangle1[0] != -1) {
                if (triangles[neighborsTriangle1[0]].neighbors[i] == triangle1Index) {
                    triangles[neighborsTriangle1[0]].neighbors[i] = triangle2Index;
                    break;
                }
            }
        }
    }

    if (adjacentSide2==0) {
        triangle1.neighbors[1] = neighborsTriangle2[1];
        triangle2.neighbors[0] = neighborsTriangle2[2];

        //aggiorniamo le adiacenze dei triangoli vicini a triangle2
        for (int i=0; i<3; i++) {
            if (neighborsTriangle2[1] != -1) {
                if (triangles[neighborsTriangle2[1]].neighbors[i] == triangle2Index) {
                    triangles[neighborsTriangle2[1]].neighbors[i] = triangle1Index;
                    break;
                }
            }
        }
    }

    else if (adjacentSide2==1) {
        triangle1.neighbors[1] = neighborsTriangle2[2];
        triangle2.neighbors[0] = neighborsTriangle2[0];

        //aggiorniamo le adiacenze dei triangoli vicini a triangle2
        for (int i=0; i<3; i++) {
            if (neighborsTriangle2[2] != -1) {
                if (triangles[neighborsTriangle2[2]].neighbors[i] == triangle2Index) {
                    triangles[neighborsTriangle2[2]].neighbors[i] = triangle1Index;
                    break;
                }
            }
        }
    }

    else if (adjacentSide2==2) {
        triangle1.neighbors[1] = neighborsTriangle2[0];
        triangle2.neighbors[0] = neighborsTriangle2[1];

        //aggiorniamo le adiacenze dei triangoli vicini a triangle2
        for (int i=0; i<3; i++) {
            if (neighborsTriangle2[0] != -1) {
                if (triangles[neighborsTriangle2[0]].neighbors[i] == triangle2Index) {
                    triangles[neighborsTriangle2[0]].neighbors[i] = triangle1Index;
                    break;
                }
            }
        }
    }

    triangle1.neighbors[2] = triangle2Index;
    triangle2.neighbors[2] = triangle1Index;

    //aggiorniamo triangles in modo vengano aggiornate le adiacenze dopo il flip
    triangles[triangle1Index] = triangle1;
    triangles[triangle2Index] = triangle2;

    //aggiornare hullTrianglesIndices solo nel caso in cui il triangolo che ha generato flip era nei triangoli del bordo
    int nHullEdgeTriangle1 = 0;  //numero di lati al bordo
    int nHullEdgeTriangle2 = 0;
    for (unsigned int k=0; k<3; k++) {
        if (triangle1.neighbors[k]==-1)
            nHullEdgeTriangle1++;
        if (triangle2.neighbors[k]==-1)
            nHullEdgeTriangle2++;
    }

    if (nHullEdgeTriangle1>=1 && nEdgeOLDTriangle1==0)           //ora questo triangolo sta nel bordo ma prima no
        hullTrianglesIndices.push_back(triangle1Index);
    else if (nHullEdgeTriangle1==0 && nEdgeOLDTriangle1>=1)      //ora questo triangolo non sta nel bordo ma prima sì
        hullTrianglesIndices.erase(remove(hullTrianglesIndices.begin(), hullTrianglesIndices.end(), triangle1Index), hullTrianglesIndices.end());

    if (nHullEdgeTriangle2>=1 && nEdgeOLDTriangle2==0)           //ora questo triangolo sta nel bordo ma prima no
        hullTrianglesIndices.push_back(triangle2Index);
    else if (nHullEdgeTriangle2==0 && nEdgeOLDTriangle2>=1)      //ora questo triangolo non sta nel bordo ma prima sì
        hullTrianglesIndices.erase(remove(hullTrianglesIndices.begin(), hullTrianglesIndices.end(), triangle2Index), hullTrianglesIndices.end());
}

 vector<Segment> Delaunay::drawSegments(vector<Triangle>& triangles)
 {
     vector<Segment> segments;
     Segment s1 = Segment(triangles[0].p1, triangles[0].p2);
     Segment s2 = Segment(triangles[0].p2, triangles[0].p3);
     Segment s3 = Segment(triangles[0].p3, triangles[0].p1);
     segments.push_back(s1);
     segments.push_back(s2);
     segments.push_back(s3);

     int trianglesSize = triangles.size();
     for (int i=1; i<trianglesSize; i++) {
         //mi creo i tre segmenti del triangolo
         Segment s1 = Segment(triangles[i].p1, triangles[i].p2);
         Segment s2 = Segment(triangles[i].p2, triangles[i].p3);
         Segment s3 = Segment(triangles[i].p3, triangles[i].p1);

         //verifico se è già presente in segments
         bool alreadyExisting1 = false;
         int drawSegmentsSize = segments.size();
         for (int k = 0; k<drawSegmentsSize; k++) {
             if (areSameSegment(s1, segments[k])) {
                     alreadyExisting1 = true;
                     break;
             }
         }
         //se non è presente, lo aggiungo
         if (alreadyExisting1 == false)
             segments.push_back(s1);

         bool alreadyExisting2 = false;
         drawSegmentsSize = segments.size();
         for (int k = 0; k<drawSegmentsSize; k++) {
             if (areSameSegment(s2, segments[k])) {
                     alreadyExisting2 = true;
                     break;
             }
         }
         if (alreadyExisting2 == false)
             segments.push_back(s2);

         bool alreadyExisting3 = false;
         drawSegmentsSize = segments.size();
         for (int k = 0; k<drawSegmentsSize; k++) {
             if (areSameSegment(s3, segments[k])) {
                     alreadyExisting3 = true;
                     break;
             }
         }
         if (alreadyExisting3 == false)
             segments.push_back(s3);
     }

     return segments;
 }

}
