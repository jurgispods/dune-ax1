FROM pederpansen/dune-ax1:dev AS build

#FROM braintwister/ubuntu-16.04
FROM debian:stable-slim

COPY --from=build /dune/dune-ax1/bin /dune/dune-ax1/bin
COPY --from=build /dune/dune-ax1/src/acme2_cyl_par /dune/dune-ax1/src/acme2_cyl_par
COPY --from=build /dune/dune-ax1/src/simulation_states /dune/dune-ax1/src/simulation_states
COPY --from=build /dune/dune-ax1/src/config /dune/dune-ax1/src/config
COPY --from=build /usr/local /usr/local
COPY --from=build /SuperLU_4.3 /SuperLU_4.3
RUN ldconfig
RUN mkdir -p /dune/output
RUN mkdir -p /dune/simulation_states

ADD bin/run.sh /dune/run.sh

ENTRYPOINT ["sh", "-c", "/dune/run.sh"]