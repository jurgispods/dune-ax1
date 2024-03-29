FROM pederpansen/dune-ax1:dev AS build

#FROM braintwister/ubuntu-16.04
FROM debian:stable-slim

COPY --from=build /dune/dune-ax1/src/acme2_cyl_par /dune/dune-ax1/src/acme2_cyl_par
COPY --from=build /dune/dune-ax1/src/acme2_cyl_par_myelin /dune/dune-ax1/src/acme2_cyl_par_myelin
COPY --from=build /usr/local /usr/local
COPY --from=build /etc/alternatives/libblas.so.3 /etc/alternatives/libblas.so.3
COPY --from=build /usr/lib/libblas.so.3 /usr/lib/libblas.so.3
COPY --from=build /usr/lib/libblas/libblas.so.3 /usr/lib/libblas//libblas.so.3
COPY --from=build /usr/lib/libblas/libblas.so.3.6.0 /usr/lib/libblas//libblas.so.3.6.0
COPY --from=build /dune/SuperLU_4.3 /dune/SuperLU_4.3
RUN ldconfig

ADD src/simulation_states /dune/dune-ax1/src/simulation_states
ADD src/config /dune/dune-ax1/src/config
RUN mkdir -p /dune/output
RUN mkdir -p /dune/simulation_states

ADD bin/run.sh /dune/run.sh

ENTRYPOINT ["/dune/run.sh"]