/** Splits a name at lowercase-to-uppercase transitions so long camelCase names can wrap. */
export const breakCamelCase = (name: string): string[] => name.split(/(?<=[a-z])(?=[A-Z])/);
