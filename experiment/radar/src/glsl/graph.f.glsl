varying vec4 graph_coord;

void main(void) {
	gl_FragColor = graph_coord / 2.0 + 0.5;
	//gl_FragColor = graph_coord / 2.0 + 0.5;
	//gl_FragColor = sin(graph_coord) / graph_coord;
}
