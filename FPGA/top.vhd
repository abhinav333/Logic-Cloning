library IEEE;
use IEEE.STD_LOGIC_1164.ALL;

entity top is
generic (N : integer; M : integer);
port(
topclk : in std_logic;
topa : in std_logic_vector(N-1 downto 0);
topb : in std_logic_vector(M-1 downto 0);
topd : in std_logic_vector(N downto 0);
topp : out std_logic_vector(N+M-1 downto 0)
);
end top;

architecture Behavioral of top is

signal reg_in_a : std_logic_vector(N-1 downto 0);
signal reg_in_b : std_logic_vector(M-1 downto 0);
signal reg_in_d: std_logic_vector(N downto 0);
signal reg_out_p : std_logic_vector(N+M-1 downto 0);

begin

--mul : entity work.optimised_accurate
--generic map(
--    N => N,
--    M => M
--)
--port map(
--    a => reg_in_a,
--    b => reg_in_b,
--    d => reg_in_d,
--    p => reg_out_p
--);

--mul_app : entity work.approximate_multiplier
--generic map(
--    N => N
-- )
--port map(
--    a => reg_in_a,
--    b => reg_in_b,
--    d => reg_in_d,
--    p => reg_out_p
--);

mul_app : entity work.mult_gen_0
port map(
    A => reg_in_a,
    B => reg_in_b,
    P => reg_out_p
);


RegProc: process(topclk) --register for the multiplier IO
begin
    if rising_edge(topclk) then
        reg_in_a <= topa;
        reg_in_b <= topb;
        reg_in_d <= topd;
        topp <= reg_out_p;
    end if;
end process RegProc;

end Behavioral;
