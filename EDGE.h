class Edge {
public:
    Edge() {}
    Edge(int a, int b, double l){
        first_vertex = a;
        second_vertex = b;
        len = l;
    }
    int first_vertex;
    int second_vertex;
    double len;
};
bool comp(const Edge& e1, const Edge& e2)
{
    return e1.len < e2.len;
}