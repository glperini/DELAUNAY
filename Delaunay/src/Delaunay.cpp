#include "Delaunay.hpp"

using namespace ProjectLibrary;


namespace ProjectLibrary{

//static constexpr double tolerance = numeric_limits<double>::epsilon();
static constexpr double tolerance = 1e-12;

//getters & setters
vector<Point> Delaunay::get_points(){
    return points;
}

Point Delaunay::get_point(int index){
    return points[index];
}

vector<Triangle> Delaunay::get_triangles(){
    return triangles;
}

Triangle Delaunay::get_triangle(int index){
    return triangles[index];
}

vector<Segment> Delaunay::get_convexHull(){
    return convexHull;
}

Segment Delaunay::get_convexHull_s(int index){
    return convexHull[index];
}

vector<Point> Delaunay::get_hullPoints(){
    return hullPoints;
}

Point Delaunay::get_hullPoint(int index){
    return hullPoints[index];
}

vector<int> Delaunay::get_hullTrianglesIndices(){
    return hullTrianglesIndices;
}

int Delaunay::get_hullTrianglesIndex(int index){
    return hullTrianglesIndices[index];
}

int Delaunay::get_triangleIndex(){
    return triangleIndex;
}

void Delaunay::set_points(vector<Point> points){
    this -> points = points;
}

void Delaunay::set_point(Point p, int index){
    this -> points[index] = p;
}

void Delaunay::set_triangles(vector<Triangle> triangles){
    this -> triangles = triangles;
}

void Delaunay::set_triangle(Triangle t, int index){
    this -> triangles[index] = t;
}

void Delaunay::set_convexHull(vector<Segment> convexHull){
    this -> convexHull = convexHull;
}

void Delaunay::set_convexHull_s(Segment s, int index){
    this -> convexHull[index] = s;
}

void Delaunay::set_hullPoints(vector<Point> hullPoints){
    this -> hullPoints = hullPoints;
}

void Delaunay::set_hullPoint(Point p, int index){
    this -> hullPoints[index] = p;
}

void Delaunay::set_hullTrianglesIndices(vector<int> hullTrianglesIndices){
    this -> hullTrianglesIndices = hullTrianglesIndices;
}

void Delaunay::set_hullTrianglesIndex(int index_s, int index){
    this -> hullTrianglesIndices[index] = index_s;
}

void Delaunay::set_triangleIndex(int triangleIndex){
    this -> triangleIndex = triangleIndex;
}


// metodi
void Delaunay::firstTriangle(vector<Point>& sortedX, vector<Point>& sortedY) {

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
    if (!p1.areSamePoint(p2) && !p1.areSamePoint(p3) && !p1.areSamePoint(p4)) {
        distinctPoints.push_back(p1);
    }
    if (!p2.areSamePoint(p3) && !p2.areSamePoint(p4)) {
        distinctPoints.push_back(p2);
    }
    if (!p3.areSamePoint(p4)) {
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
        tryTriangles.push_back(Triangle({p1, p2, p3}));
        tryTriangles.push_back(Triangle({p1, p2, p4}));
        tryTriangles.push_back(Triangle({p2, p3, p4}));
        tryTriangles.push_back(Triangle({p1, p3, p4}));

        //calcolo le aree dei triangoli e mi tengo la piu grande
        double bestArea = 0;
        double areaTriangle;

        for (unsigned int i=0; i<4; i++) {
            areaTriangle = tryTriangles[i].Area();
            if (areaTriangle>bestArea){
                bestArea = areaTriangle;
                bestTriangle = tryTriangles[i];
            }
        }
    }

    //se ho 3 punti distinti, basta costruire il triangolo grande con quei 3 punti
    else if (distinctPoints.size() == 3) {
        bestTriangle = Triangle({distinctPoints[0], distinctPoints[1], distinctPoints[2]});
    }

    //se ho 2 punti distinti:
    else if (distinctPoints.size() == 2) {
        Point newp3 = sortedX[1];  //prendo il secondo punto più piccolo secondo x
        Point newp4 = sortedX[sortedX.size() - 2];  //prendo il penultimo punto più piccolo secondo x
        Triangle triangle3 = Triangle({distinctPoints[0], distinctPoints[1], newp3});
        Triangle triangle4 = Triangle({distinctPoints[0], distinctPoints[1], newp4});
        double areap3 = triangle3.Area();  //area con il secondo punto
        double areap4 = triangle4.Area();  //area con il penultimo punto
        if (areap3<areap4)
            bestTriangle = triangle4;
        else
            bestTriangle = triangle3;
    }

    //tolgo da sortedX i 3 vertici del triangolo grande appena creato
    int p1_index = searchPoint(sortedX, bestTriangle.get_point(0));
    int p2_index = searchPoint(sortedX, bestTriangle.get_point(1));
    int p3_index = searchPoint(sortedX, bestTriangle.get_point(2));

    if (p1_index!=-1 && p2_index!=-1 && p3_index!=-1) {
        sortedX.erase(sortedX.begin() +  bestTriangle.get_point(0).get_id());
        if (p1_index<p2_index) {
            sortedX.erase(sortedX.begin() +  bestTriangle.get_point(1).get_id() - 1);
            if (p2_index<p3_index)
                sortedX.erase(sortedX.begin() +  bestTriangle.get_point(2).get_id() - 2);
            else
                sortedX.erase(sortedX.begin() +  bestTriangle.get_point(2).get_id() - 1);
        }
        else { //(p1_index>p2_index)
            sortedX.erase(sortedX.begin() +  bestTriangle.get_point(1).get_id());
            if (p1_index<p3_index)
                sortedX.erase(sortedX.begin() +  bestTriangle.get_point(2).get_id() - 1);
            else //(p1_index>p3_index)
                sortedX.erase(sortedX.begin() +  bestTriangle.get_point(2).get_id());
        }
    }
    else
        cerr << "Errore" << endl;

    //aggiungo questo triangolo ordinato in senso antiorario al vettore dei triangoli
    bestTriangle.antiClockWiseOrder();
    triangles.push_back(bestTriangle);

    //agigungo i punti al guscio
    hullPoints.push_back(bestTriangle.get_point(0));
    hullPoints.push_back(bestTriangle.get_point(1));
    hullPoints.push_back(bestTriangle.get_point(2));

    //aggiungo segmenti di bordo
    convexHull.push_back(Segment({bestTriangle.get_point(0), bestTriangle.get_point(1)}));
    convexHull[0].lineCoefficients();

    convexHull.push_back(Segment({bestTriangle.get_point(1), bestTriangle.get_point(2)}));
    convexHull[1].lineCoefficients();

    convexHull.push_back(Segment({bestTriangle.get_point(2), bestTriangle.get_point(0)}));
    convexHull[2].lineCoefficients();

    //neighborhood settato di default
    hullTrianglesIndices.push_back(0);
}

//Verifica se point è interno/esterno al circumcerchio di ogni traingolo.
//Si ferma o quando trova un circumcerchio a cui è interno o quando li ha verificati tutti e sta fuori
int Delaunay::isPointInsideCircumcircle(Point point) {

    int triangleSize = triangles.size();
    for (int i = triangleIndex; i<triangleSize; i++) {
        //seleziona triangolo
        Triangle triangle = triangles[i];


        double ax = triangle.get_point(0).get_x() - point.get_x(); //Ax-Dx
        double ay = triangle.get_point(0).get_y() - point.get_y(); //Ay-Dy
        double bx = triangle.get_point(1).get_x() - point.get_x(); //Bx-Dx
        double by = triangle.get_point(1).get_y() - point.get_y(); //By-Dy
        double cx = triangle.get_point(2).get_x() - point.get_x(); //Cx-Dx
        double cy = triangle.get_point(2).get_y() - point.get_y(); //Cy-Dy

        double d = ax * (by * (pow(cx,2) + pow(cy,2)) - cy * (pow(bx,2) + pow(by,2))) -
                bx * (ay * (pow(cx,2) + pow(cy,2)) - cy * (pow(ax,2) + pow(ay,2))) +
                cx * (ay * (pow(bx,2) + pow(by,2)) - by * (pow(ax,2) + pow(ay,2)));

        if (d > 0)   //d>0 dice che il punto sta dentro circumucerchio del triangolo
            return i;  //restituisce indice di triangolo appena testato
    }
    return triangleIndex = -1;
}

//Verifica se un punto è interno/esterno al triangolo rilevato dal circumcerchio di isPointInsideCircumcircle
int Delaunay::isPointInsideTriangle(Point point) {

    //prendo il triangolo dal vettore triangles
    Triangle triangle = triangles[triangleIndex];

    //determinante che serve per determinare se il punto è interno al triangolo
    double d1 = (point.get_x() - triangle.get_point(0).get_x()) * (triangle.get_point(2).get_y() - triangle.get_point(0).get_y()) - (point.get_y() - triangle.get_point(0).get_y()) * (triangle.get_point(2).get_x() - triangle.get_point(0).get_x());
    double d2 = (triangle.get_point(0).get_x() - point.get_x()) * (triangle.get_point(1).get_y() - triangle.get_point(0).get_y()) - (triangle.get_point(0).get_y() - point.get_y()) * (triangle.get_point(1).get_x() - triangle.get_point(0).get_x());
    double d3 = (triangle.get_point(1).get_x() - triangle.get_point(0).get_x()) * (triangle.get_point(2).get_y() - triangle.get_point(0).get_y()) - (triangle.get_point(1).get_y() - triangle.get_point(0).get_y()) * (triangle.get_point(2).get_x() - triangle.get_point(0).get_x());
    double a = d1/d3;
    double b = d2/d3;

    //if (a>0 && b>0 && (a+b)<1)    //se tutti 3 determinanti sono positivi, si trova all’interno del triangolo
    if (a > 0 && b > 0 && (a + b) < 1)
        return triangleIndex;       //restituisce triangleIndex (che davamo come input alla funzione) se il punto è interno al triangolo
    else
        return -1;        //altrimenti -1 se è esterno alla triangolazione
}

int Delaunay::findTriangleContainingPoint(Point point) {

    if (isPointInsideCircumcircle(point)!=-1) {  //se il punto è interno al circumcerchio
        if (isPointInsideTriangle(point)!=-1)  //se il punto è interno al triangolo
            return triangleIndex;
        else {  //se il punto è esterno al triangolo ripeto
            int triangleSizeMinus1 = triangles.size()-1;
            if (triangleIndex<triangleSizeMinus1) {
                triangleIndex = triangleIndex + 1;
                findTriangleContainingPoint(point);
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
void Delaunay::inTriangle(Point point){

    //prendo il triangolo dal vettore triangles
    Triangle triangle = triangles[triangleIndex];

    //salvo i punti del triangolo nelle variabili v1, v2, v3
    Point v1 = triangle.get_point(0);
    Point v2 = triangle.get_point(1);
    Point v3 = triangle.get_point(2);

    //creo 3 nuovi triangoli vuoti
    Triangle t1, t2, t3;
    int triangleIndex2 = triangles.size();
    int triangleIndex3 = triangles.size() + 1;

    t1.set_points({v1, v2, point}); //primo vertice del triangolo grande, secondo vertice del triangolo grande, terzo vertice è il punto interno al triangolo grande
    t1.set_neighbors({triangle.get_neighbor(0), triangleIndex2, triangleIndex3});   //primo lato + l'indice del triangolo vicino del triangolo grande che ha in comune il lato con t1
    //secondo lato + indice di t2, che subito dopo creerò
    //terzo lato + indice di t3, che subito dopo creerò

    //le adiacenze del "vecchio" triangolo adiacente a t1 rimangono invariate

    t2.set_points({v2, v3, point}); //secondo vertice del triangolo grande, terzo vertice del triangolo grande, il terzo vertice è il punto interno al triangolo grande
    t2.set_neighbors({triangle.get_neighbor(1), triangleIndex3, triangleIndex});    //primo lato + l'indice del triangolo vicino del triangolo grande che ha in comune il lato con t2
    //secondo lato + indice di t3
    //terzo lato + indice di t1

    //aggiorniamo le adiacenze del triangolo vicino a t2
    for (int i=0; i<3; i++) {
        if (t2.get_neighbor(0) != -1) {
            if (triangles[t2.get_neighbor(0)].get_neighbor(i) == triangleIndex) {
                triangles[t2.get_neighbor(0)].set_neighbor(triangleIndex2, i);
                break;
            }
        }
    }

    t3.set_points({v3, v1, point});
    t3.set_neighbors({triangle.get_neighbor(2), triangleIndex, triangleIndex2});    //primo lato + indice del triangolo vicino del triangolo grande che ha in comune il lato con t3
    //secondo lato + indice di t1
    //terzo lato + indice di t2

    //aggiorniamo le adiacenze del triangolo vicino a t3
    for (int i=0; i<3; i++) {
        if (t3.get_neighbor(0) != -1) {
            if (triangles[t3.get_neighbor(0)].get_neighbor(i) == triangleIndex) {
                triangles[t3.get_neighbor(0)].set_neighbor(triangleIndex3, i);
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
    if (t1.get_neighbor(0)==-1 || t2.get_neighbor(0)==-1 || t3.get_neighbor(0)==-1) {
        hullTrianglesIndices.erase(remove(hullTrianglesIndices.begin(), hullTrianglesIndices.end(), triangleIndex), hullTrianglesIndices.end());
        if (t1.get_neighbor(0)==-1)  //se t1 non ha vicini, inserisco il suo indice (che adesso è triangleIndex) dentro hullTrianglesIndices
            hullTrianglesIndices.push_back(triangleIndex);
        if (t2.get_neighbor(0)==-1)  //se t2 non ha vicini, inserisco il suo indice (che corrisponde al penultimo elemento di triangles) dentro hullTrianglesIndices
            hullTrianglesIndices.push_back(triangleIndex2);
        if (t3.get_neighbor(0)==-1)  //se t3 non ha vicini, inserisco il suo indice (che corrisponde all'ultimo elemento di triangles) dentro hullTrianglesIndices
            hullTrianglesIndices.push_back(triangleIndex3);
    }

    verifyDelaunayCondition(t1);
    triangleIndex = triangleIndex2;
    verifyDelaunayCondition(t2);
    triangleIndex = triangleIndex3;
    verifyDelaunayCondition(t3);

}

//Caso punto esterno triangolazione
//unisce punto esterno con tutti i vertici su guscio e rimuove segmenti che hanno intersezione
//aggiorna neighbors e guscio
void Delaunay::outTriangle(Point outsidePoint){

    //creo un vettore che conterrà i nuovi segmenti
    vector<Segment> newSegments;
    Segment temporary = Segment();  //creazione di segmento temporaneo dove salvo i coefficienti
    for (unsigned int i=0; i<hullPoints.size(); i++){  //cicliamo su tutti i punti del guscio
        //creiamo una retta (vettore con m e q) dato il punto esterno e l'i-esimo punto del guscio
        temporary.set_points({outsidePoint, hullPoints[i]});
        temporary.lineCoefficients();
        Vector2d possibleNewLine = temporary.get_coefficients();
        /*iniziamo a contare il numero di intersezioni tra l'i-esima retta (corrispontende all'i-esimo punto)
            e le rette che corrispondono ai segmenti del guscio*/
        int numberOfIntersections = 0;
        for (unsigned int j=0; j<convexHull.size(); j++){
            //troviamo il punto di intersezione tra l'i-esima retta e la j-esima retta (corrispondente al j-esimo segmento del guscio):
            Point intersection = findIntersection(possibleNewLine, convexHull[j].get_coefficients());
            //se almeno una delle due coordinate esiste:
            if (!isnan(intersection.get_x()) || !isnan(intersection.get_y())){
                //calcolo la lunghezza del j-esimo segmento del guscio
                double lenghtSegment = convexHull[j].get_point(0).getDistance(convexHull[j].get_point(1));
                /*se entrambe le distanze tra il punto intersezione e i due punti del j-esimo segmento
                    sono minori della lunghezza del j-esimo segmento, ho intersezione e devo aggiungerla a numberOfIntersections:*/

                double distance1 = intersection.getDistance(convexHull[j].get_point(0));
                double distance2 = intersection.getDistance(convexHull[j].get_point(1));

                //se m e q delle due rette confrontate per intersezione sono uguali, allora le rette sono coincidenti. l'intersezione va ritenuta valida
                if (areEqual(possibleNewLine, convexHull[j].get_coefficients()))
                    numberOfIntersections++;
                else if (distance1<=(lenghtSegment+tolerance) && distance2<=(lenghtSegment+tolerance)) {
                    double lengthSide = outsidePoint.getDistance(hullPoints[i]);
                    double distance = outsidePoint.getDistance(intersection);
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
        if (numberOfIntersections==2){
            Segment temporarySegment = Segment({outsidePoint, hullPoints[i]});
            temporarySegment.set_coefficients(possibleNewLine);
            newSegments.push_back(temporarySegment);
        }
    }

    //aggiungo il nuovo punto a hullPoints
    hullPoints.push_back(outsidePoint);

    int numberOfNewTriangles = 0;
    for (unsigned int i=0; i<newSegments.size(); i++){
        for (unsigned int j=0; j<newSegments.size(); j++){
            if (i<j) {
                //creo un segmento con due punti della convexHull
                Segment verifySegment = Segment({newSegments[i].get_point(1), newSegments[j].get_point(1)});

                bool existed = false;
                for (unsigned int q=0; q<convexHull.size(); q++) {
                    if (verifySegment.areSameSegment(convexHull[q])) {
                        existed = true;
                        break;
                    }
                }

                //verifico se segmento che chiude triangolo costruito con i due segmenti che sto accoppiando fa parte del guscio, qui == quindi SI
                if (existed) {
                    Triangle newTriangle = Triangle({outsidePoint, newSegments[i].get_point(1), newSegments[j].get_point(1)});
                    newTriangle.antiClockWiseOrder();
                    int newTriangleIndex = triangles.size();
                    triangles.push_back(newTriangle);
                    numberOfNewTriangles = numberOfNewTriangles+1;

                    //aggiorno le adiacenze
                    unsigned int temporary = 0;
                    Triangle hullTriangle;
                    int hullTriangleIndicesSize = hullTrianglesIndices.size();
                    for (int k=0; k<hullTriangleIndicesSize; k++){                                                     //per ogni triangolo del guscio controllo se il segmento di adiacenza con il nuovo triangolo aggiunto gli appartiene (newTriangle.get_point(0) = outsidePoint)
                        hullTriangle = triangles[hullTrianglesIndices[k]];
                        if (newTriangle.get_point(1).areSamePoint(hullTriangle.get_point(0)) &&                                        //confronto questi vertici di hullTriangle perchè so che sono ordinati in senso antiorario
                                newTriangle.get_point(2).areSamePoint(hullTriangle.get_point(2)))  {
                            hullTriangle.set_neighbor(newTriangleIndex, 2);                                    //aggiorno l'indice del vicino mettendogli l'indice del triangolo appena aggiunto
                            temporary = k;
                            break;
                        }

                        else if (newTriangle.get_point(1).areSamePoint(hullTriangle.get_point(1)) &&                                        //confronto questi vertici di hullTriangle perchè so che sono ordinati in senso antiorario
                                 newTriangle.get_point(2).areSamePoint(hullTriangle.get_point(0))) {
                            hullTriangle.set_neighbor(newTriangleIndex, 0);                                    //aggiorno l'indice del vicino mettendogli l'indice del triangolo appena aggiunto
                            temporary = k;
                            break;
                        }
                        else if (newTriangle.get_point(1).areSamePoint(hullTriangle.get_point(2)) &&                                       //confronto questi vertici di hullTriangle perchè so che sono ordinati in senso antiorario
                                 newTriangle.get_point(2).areSamePoint(hullTriangle.get_point(1))) {
                            hullTriangle.set_neighbor(newTriangleIndex, 1);                                    //aggiorno l'indice del vicino mettendogli l'indice del triangolo appena aggiunto
                            temporary = k;
                            break;
                        }
                        //se non era indice k giusto, esco e continuo con il successivo
                    }
                    triangles[hullTrianglesIndices[temporary]].set_neighbors(hullTriangle.get_neighbors());
                    newTriangle.set_neighbor(hullTrianglesIndices[temporary], 1);                                       //ora aggiorno il vicino del nuovo triangolo dal lato adiacente con un vecchio
                    triangles[newTriangleIndex].set_neighbors(newTriangle.get_neighbors());

                    //controlliamo se il triangolo adiacente a quello creato ha lati sul bordo
                    bool onSide = false;
                    for (int l = 0; l<3; l++) {
                        if (triangles[hullTrianglesIndices[temporary]].get_neighbor(l) == -1) {
                            onSide = true;
                            break;
                        }
                    }
                    //se non ha più lati sul bordo, eliminiamo l'indice da hullTrianglesIndices
                    if (onSide == false)
                        hullTrianglesIndices.erase(remove(hullTrianglesIndices.begin(), hullTrianglesIndices.end(), hullTrianglesIndices[temporary]), hullTrianglesIndices.end());      //rimuovo l'indice del triangolo adiacente che prima era di bordo e ora non è più

                    //elimino il lato che prima era nella convexHull e ora non lo è più
                    Segment oldHullSegment = Segment({newTriangle.get_point(1), newTriangle.get_point(2)});
                    int convexHullSize = convexHull.size();
                    for (int k=0; k<convexHullSize; k++) {
                        if (oldHullSegment.areSameSegment(convexHull[k])) {
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
        Point p1 = triangles[trianglesSizeMinus1].get_point(0);
        Point p2 = triangles[trianglesSizeMinus1].get_point(1);
        Point p3 = triangles[trianglesSizeMinus1].get_point(2);
        Segment s1 = Segment({p1, p2});
        Segment s2 = Segment({p3, p1});
        s1.lineCoefficients();
        s2.lineCoefficients();
        convexHull.push_back(s1);
        convexHull.push_back(s2);
    }

    else {
        for (int i=0; i<numberOfNewTriangles; i++) {
            for (int j=0; j<numberOfNewTriangles; j++) {
                if (i!=j) {                                                                                     //così non confronto con se stesso lo stesso triangolo!
                    Triangle triangleToCompare1 = triangles[triangles.size()-numberOfNewTriangles + i];         //salvo temporaneamente per facilitare la lettura del codice
                    Triangle triangleToCompare2 = triangles[triangles.size()-numberOfNewTriangles + j];
                    if (triangleToCompare1.get_point(1).areSamePoint(triangleToCompare2.get_point(2))) {
                        triangleToCompare1.set_neighbor(triangles.size()-numberOfNewTriangles + j, 0);            //aggiorno in triangles l'indice del triangolo corrispondente a triangleToCompare1 e gli assegno l'indice di quello di fianco a destra (corrispondente a triangleToCompare2)
                        triangleToCompare2.set_neighbor(triangles.size()-numberOfNewTriangles + i, 2);

                        int hullPointsSize = hullPoints.size();
                        for (int k=0; k<hullPointsSize; k++) {
                            //elimino il punto in comune fra i due triangoli creati perchè non è più di bordo
                            if (triangleToCompare1.get_point(1).areSamePoint(hullPoints[k])) {
                                hullPoints.erase(next(hullPoints.begin(), k));
                                break;
                            }
                        }
                        Segment convexHullSegment;
                        if (triangleToCompare1.get_neighbor(2)==-1) {                                                //aggiungo segmento di bordo a convexHull
                            convexHullSegment.set_points({triangleToCompare1.get_point(2), outsidePoint});
                            convexHullSegment.lineCoefficients();
                            convexHull.push_back(convexHullSegment);
                        }
                        if (triangleToCompare1.get_neighbor(0)==-1) {
                            convexHullSegment.set_points({outsidePoint, triangleToCompare1.get_point(1)});
                            convexHullSegment.lineCoefficients();
                            convexHull.push_back(convexHullSegment);
                        }

                        if (triangleToCompare2.get_neighbor(2)==-1) {                                                //aggiungo segmento di bordo a convexHull
                            convexHullSegment.set_points({triangleToCompare2.get_point(2), outsidePoint});
                            convexHullSegment.lineCoefficients();
                            convexHull.push_back(convexHullSegment);
                        }
                        if (triangleToCompare2.get_neighbor(0)==-1) {
                            convexHullSegment.set_points({outsidePoint, triangleToCompare2.get_point(1)});
                            convexHullSegment.lineCoefficients();
                            convexHull.push_back(convexHullSegment);
                        }

                        triangles[triangles.size()-numberOfNewTriangles + i].set_neighbor(triangleToCompare1.get_neighbor(0), 0);
                        triangles[triangles.size()-numberOfNewTriangles + j].set_neighbor(triangleToCompare2.get_neighbor(2), 2);
                    }
                }
            }
        }
    }

    int firstNewTriangleIndex = triangles.size()-numberOfNewTriangles;
    for (int i=0; i<numberOfNewTriangles; i++) {
        int triangleIndex = firstNewTriangleIndex+i;
        for (int j=0; j<3; j++) {
            if (triangles[triangleIndex].get_neighbor(j) == -1) {
                hullTrianglesIndices.push_back(triangleIndex);
                break;
            }
        }
    }

    for (int i=0; i<numberOfNewTriangles; i++) {
        triangleIndex = firstNewTriangleIndex+i;
        verifyDelaunayCondition(triangles[triangleIndex]);
    }
}

//Verifico la condizione di Delunay
void Delaunay::verifyDelaunayCondition(Triangle tr){

    vector<int> flippedTrianglesIndices;
    Triangle triangle = tr;
    double angle1;
    double angle2;
    Triangle nearTriangle;
    Segment temporarySegment;
    //i sarà l'indice che identifica il lato del triangolo 1 su cui ho verificato avere adiacenza con triangolo 2
    //j sarà l'indice che identifica il lato del triangolo 2 su cui ho verificato avere adiacenza con triangolo 1
    for (int i=0; i<3; i++) {                    //per ogni lato
        if (triangle.get_neighbor(i)!=-1) {      //triangolo testato ha un triangolo adiacente su quel lato
            nearTriangle = triangles[triangle.get_neighbor(i)];             //il triangolo vicino è quello indicato dall'indice contenuto in neighbors
            for (int j=0; j<3; j++) {
                if (nearTriangle.get_neighbor(j) == triangleIndex) {         //capisco su quale lato del triangolo vicino i due triangoli sono adiacenti
                    if (i==0) {
                        Segment segment1 = Segment({triangle.get_point(0), triangle.get_point(2)});
                        Segment segment2 = Segment({triangle.get_point(1), triangle.get_point(2)});
                        angle1 = segment1.getAngle(segment2);
                    }
                    //angle1 = getAngle(Segment(triangle.get_point(0), triangle.get_point(2)), Segment(triangle.get_point(1), triangle.get_point(2)));                   //calcolo l'angolo opposto al lato adiacente nel triangolo di partenza
                    else if (i==1) {
                        temporarySegment = Segment({triangle.get_point(1), triangle.get_point(0)});
                        angle1 = temporarySegment.getAngle(Segment({triangle.get_point(2), triangle.get_point(0)}));
                    }
                    else if (i==2) {
                        temporarySegment = Segment({triangle.get_point(2), triangle.get_point(1)});
                        angle1 = temporarySegment.getAngle(Segment({triangle.get_point(0), triangle.get_point(1)}));
                    }
                    if (j==0) {
                        temporarySegment = Segment({nearTriangle.get_point(0), nearTriangle.get_point(2)});
                        angle2 = temporarySegment.getAngle(Segment({nearTriangle.get_point(1), nearTriangle.get_point(2)}));         //calcolo l'angolo opposto al lato adiacente nel triangolo vicino
                    }
                    else if (j==1) {
                        temporarySegment = Segment({nearTriangle.get_point(1), nearTriangle.get_point(0)});
                        angle2 = temporarySegment.getAngle(Segment({nearTriangle.get_point(2), nearTriangle.get_point(0)}));
                    }
                    else if (j==2) {
                        temporarySegment = Segment({nearTriangle.get_point(2), nearTriangle.get_point(1)});
                        angle2 = temporarySegment.getAngle(Segment({nearTriangle.get_point(0), nearTriangle.get_point(1)}));
                    }
                    //se somma è maggiore di 180° faccio flip
                    if (angle1 + angle2 > 180.0) {
                        int triangle2Index = triangle.get_neighbor(i);
                        flippedTrianglesIndices.push_back(triangle2Index);
                        flipTriangles(triangleIndex, triangle2Index, i, j);
                        triangle = triangles[triangleIndex];
                        i = 0;
                        break;
                    }

                }
            }
        }
    }
    int numberOfFlips = flippedTrianglesIndices.size();
    for (int k=0; k<numberOfFlips; k++) {
        triangleIndex = flippedTrianglesIndices[k];
        verifyDelaunayCondition(triangles[triangleIndex]);
    }
}

//Flippo sul triangolo
//l'ordine con cui passo gli indici dei triangoli da flippare è ininfluente
void Delaunay::flipTriangles(int triangle1Index, int triangle2Index, int adjacentSide1, int adjacentSide2){

    Triangle oldTriangle1 = triangles[triangle1Index];
    Triangle oldTriangle2 = triangles[triangle2Index];
    //salvo quanti lati sul bordo avevano triangle1 e triangle2 prima di essere flippati. Mi servirà per dopo
    int nEdgeOLDTriangle1 = 0;
    int nEdgeOLDTriangle2 = 0;
    for (int k=0; k<3; k++) {
        if (oldTriangle1.get_neighbor(k)==-1)
            nEdgeOLDTriangle1++;
        if (oldTriangle2.get_neighbor(k)==-1)
            nEdgeOLDTriangle2++;
    }

    //in base al lato che ha generato adiacenza salvo i valori dei punti
    Point v1Common;
    Point v2Common;
    Point v3FromTriangle1;
    Point v3FromTriangle2;

    if (adjacentSide1==0) {
        v1Common = oldTriangle1.get_point(0);
        v2Common = oldTriangle1.get_point(1);
        v3FromTriangle1 = oldTriangle1.get_point(2);
    }
    else if (adjacentSide1==1) {
        v1Common = oldTriangle1.get_point(1);
        v2Common = oldTriangle1.get_point(2);
        v3FromTriangle1 = oldTriangle1.get_point(0);
    }
    else if (adjacentSide1==2) {
        v1Common = oldTriangle1.get_point(2);
        v2Common = oldTriangle1.get_point(0);
        v3FromTriangle1 = oldTriangle1.get_point(1);
    }

    if (adjacentSide2==0)
        v3FromTriangle2 = oldTriangle2.get_point(2);
    else if (adjacentSide2==1)
        v3FromTriangle2 = oldTriangle2.get_point(0);
    else if (adjacentSide2==2)
        v3FromTriangle2 = oldTriangle2.get_point(1);

    //devo salvare prima i vicini, poi sostituisco triangolo e poi ripristino vicini
    vector<int> neighborsTriangle1 = oldTriangle1.get_neighbors();
    vector<int> neighborsTriangle2 = oldTriangle2.get_neighbors();
    Triangle triangle1 = Triangle({v3FromTriangle1, v1Common, v3FromTriangle2});       //dovrebbero già essere ordinati in senso antiorario
    Triangle triangle2 = Triangle({v3FromTriangle2, v2Common, v3FromTriangle1});

    if (adjacentSide1==0) {
        triangle1.set_neighbor(neighborsTriangle1[2], 0);
        triangle2.set_neighbor(neighborsTriangle1[1], 1);

        //aggiorniamo le adiacenze dei triangoli vicini a triangle1
        for (int i=0; i<3; i++) {
            if (neighborsTriangle1[1] != -1) {
                if (triangles[neighborsTriangle1[1]].get_neighbor(i) == triangle1Index) {
                    triangles[neighborsTriangle1[1]].set_neighbor(triangle2Index, i);
                    break;
                }
            }
        }
    }

    else if (adjacentSide1==1) {
        triangle1.set_neighbor(neighborsTriangle1[0], 0);
        triangle2.set_neighbor(neighborsTriangle1[2], 1);

        //aggiorniamo le adiacenze dei triangoli vicini a triangle1
        for (int i=0; i<3; i++) {
            if (neighborsTriangle1[2] != -1) {
                if (triangles[neighborsTriangle1[2]].get_neighbor(i) == triangle1Index) {
                    triangles[neighborsTriangle1[2]].set_neighbor(triangle2Index, i);
                    break;
                }
            }
        }
    }

    else if (adjacentSide1==2) {
        triangle1.set_neighbor(neighborsTriangle1[1], 0);
        triangle2.set_neighbor(neighborsTriangle1[0], 1);

        //aggiorniamo le adiacenze dei triangoli vicini a triangle1
        for (int i=0; i<3; i++) {
            if (neighborsTriangle1[0] != -1) {
                if (triangles[neighborsTriangle1[0]].get_neighbor(i) == triangle1Index) {
                    triangles[neighborsTriangle1[0]].set_neighbor(triangle2Index, i);
                    break;
                }
            }
        }
    }

    if (adjacentSide2==0) {
        triangle1.set_neighbor(neighborsTriangle2[1], 1);
        triangle2.set_neighbor(neighborsTriangle2[2], 0);

        //aggiorniamo le adiacenze dei triangoli vicini a triangle2
        for (int i=0; i<3; i++) {
            if (neighborsTriangle2[1] != -1) {
                if (triangles[neighborsTriangle2[1]].get_neighbor(i) == triangle2Index) {
                    triangles[neighborsTriangle2[1]].set_neighbor(triangle1Index, i);
                    break;
                }
            }
        }
    }

    else if (adjacentSide2==1) {
        triangle1.set_neighbor(neighborsTriangle2[2], 1);
        triangle2.set_neighbor(neighborsTriangle2[0], 0);

        //aggiorniamo le adiacenze dei triangoli vicini a triangle2
        for (int i=0; i<3; i++) {
            if (neighborsTriangle2[2] != -1) {
                if (triangles[neighborsTriangle2[2]].get_neighbor(i) == triangle2Index) {
                    triangles[neighborsTriangle2[2]].set_neighbor(triangle1Index, i);
                    break;
                }
            }
        }
    }

    else if (adjacentSide2==2) {
        triangle1.set_neighbor(neighborsTriangle2[0], 1);
        triangle2.set_neighbor(neighborsTriangle2[1], 0);

        //aggiorniamo le adiacenze dei triangoli vicini a triangle2
        for (int i=0; i<3; i++) {
            if (neighborsTriangle2[0] != -1) {
                if (triangles[neighborsTriangle2[0]].get_neighbor(i) == triangle2Index) {
                    triangles[neighborsTriangle2[0]].set_neighbor(triangle1Index, i);
                    break;
                }
            }
        }
    }

    triangle1.set_neighbor(triangle2Index, 2);
    triangle2.set_neighbor(triangle1Index, 2);

    //aggiorniamo triangles in modo vengano aggiornate le adiacenze dopo il flip
    triangles[triangle1Index] = triangle1;
    triangles[triangle2Index] = triangle2;

    //aggiornare hullTrianglesIndices solo nel caso in cui il triangolo che ha generato flip era nei triangoli del bordo
    int nHullEdgeTriangle1 = 0;  //numero di lati al bordo
    int nHullEdgeTriangle2 = 0;
    for (unsigned int k=0; k<3; k++) {
        if (triangle1.get_neighbor(k)==-1)
            nHullEdgeTriangle1++;
        if (triangle2.get_neighbor(k)==-1)
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

int Delaunay::searchPoint(vector<Point> points, Point target){

    int sx = 0;
    int dx = points.size();
    bool found = false;
    int index = -1;
    while ((sx<=dx) && (!found)){
        int cx = (sx+dx)/2;
        if (points[cx].get_x()>target.get_x())
            dx = cx - 1;
        else if (points[cx].get_x()<target.get_x())
            sx = cx + 1;
        else {
            if (points[cx].get_y()==target.get_y()){
                index = cx;
                found = true;
            }
            else if (points[cx].get_y()>target.get_y()){
                cx--;
                while (!found){
                    if (points[cx].get_y()==target.get_y()){
                        index = cx;
                        found = true;
                    }
                    cx--;
                }
            }
            else {
                cx++;
                while (!found){
                    if (points[cx].get_y()==target.get_y()){
                        index = cx;
                        found = true;
                    }
                    cx++;
                }
            }
        }
    }
    return index;
}

Point Delaunay::findIntersection(Vector2d r1, Vector2d r2){
    Point intersectionPoint;
    if (!isnan(r1[0]) && !isnan(r2[0])){  //se le rette non sono verticali
        if (r1[0]==r2[0]){
            if (r1[1]==r2[1])        //caso di due rette coincidenti
                intersectionPoint = Point({10, r1[0]*10 + r1[1]});     //assegnato valore casuale solo per avere punto valido
            else
                intersectionPoint = Point({numeric_limits<double>::quiet_NaN(), numeric_limits<double>::quiet_NaN()});
        }
        else {      // (r1[0]!=r2[0])
            double x = (r2[1]-r1[1])/(r1[0]-r2[0]);     //si trova la x con x = (q2-q1)/(m1-m2)
            intersectionPoint = Point({x, r1[0]*x + r1[1]});     //si sostituisce in una delle due rette per trovare la y
        }
    }
    //se invece una retta è verticale o entrambe sono verticali:
    else if (isnan(r1[0]) && isnan(r2[0]))       //se entrambe le rette sono verticali non c'è intersezione
        intersectionPoint = Point({numeric_limits<double>::quiet_NaN(), numeric_limits<double>::quiet_NaN()});
    else if (isnan(r1[0]))     //se la prima retta è verticale
        intersectionPoint = Point({r1[1], r2[0]*r1[1] + r2[1]});
    else if (isnan(r2[0]))     //se la seconda retta è verticale
        intersectionPoint = Point({r2[1], r1[0]*r2[1] + r1[1]});
    return intersectionPoint;
}

template<typename T>
bool Delaunay::areEqual(T v1, T v2) {
    if (isnan(v1[0]) && isnan(v2[0]) && v1[1] == v2[1])
        return true;
    else if (v1[0] == v2[0] && v1[1] == v2[1])
        return true;
    else
        return false;
}

}
