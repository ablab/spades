<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
    <xsl:output encoding="UTF-8" method="html" omit-xml-declaration="yes" indent="yes"/>

    <xsl:template match="tree">
        <html>
            <title>MGRA tree</title>
            <body>
                <h2>MGRA tree</h2>
                <table border="1" cellpadding="10">
                    <xsl:apply-templates select="row"/>
                </table>
            </body>
        </html>
    </xsl:template>

    <xsl:template match="row">
        <tr>
             <xsl:apply-templates select="cell"/>
        </tr>
    </xsl:template>

    <xsl:template match="cell">
        <td rowspan="{height}" colspan="{width}" align="center">
            <xsl:apply-templates select="length"/>
            <strong><xsl:value-of select="text"/></strong>
        </td>
    </xsl:template>

    <xsl:template match="length">
        <xsl:value-of select="."/>
        <br/>
    </xsl:template>
</xsl:stylesheet>
